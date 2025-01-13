#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

//Important - many of the functions/datatypes have amrex:: infront.
using namespace amrex;

static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);
    //declare ParmParse
    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(false);
    //find filename
    std::string plotFileName;
    pp.get("infile",plotFileName);

    //count number of progress variables, declare and fill array
    int numProgressVars = pp.countval("progressNames");
    Vector<std::string> progressNames(numProgressVars);
    pp.getarr("progressNames",progressNames);

    //get the burnt and unburnt values
    Vector<Real> host_c_u(numProgressVars),host_c_b(numProgressVars); 
    pp.getarr("burntProgress",host_c_b);
    pp.getarr("unburntProgress",host_c_u);

    //put it onto GPU array
    Gpu::DeviceVector<Real> device_c_u(host_c_u.size()),device_c_b(host_c_b.size());
    Gpu::copyAsync(amrex::Gpu::hostToDevice, 
                          host_c_u.begin(), 
                          host_c_u.end(), 
                          device_c_u.begin());
    Gpu::copyAsync(amrex::Gpu::hostToDevice, 
                          host_c_b.begin(), 
                          host_c_b.end(), 
                          device_c_b.begin());

    //Synchronize streams on the GPU
    Gpu::streamSynchronize();

    //Get the pointer to data in GPU memory
    Real* d_c_u = device_c_u.data();
    Real* d_c_b = device_c_b.data();
    
    //check if we're periodic for geometry information
    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }

    //initialise DataServices (IO)
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();
    //How many levels of data do we want? Default to all
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    //make a RealBox (physical geometry, rather than just cells), just for geometry information
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));
    
    //Get the in names and create outnames (destfillcomps maps innames from whatever index they are in file to 0->numPV)
    const int nCompIn = numProgressVars;
    const int nCompOut = numProgressVars;
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);

    for (int n = 0; n < nCompIn; n++) {
      inNames[n] = progressNames[n];
      destFillComps[n] = n;
    }
    for (int n = 0; n < nCompOut; n++) {
      outNames[n] = "c("+progressNames[n]+")";
    }
    //declare Vector of MultiFab, one MultiFab for each level
    Vector<MultiFab> outdata(Nlev);
    //also need to hold the geometry of each level
    Vector<Geometry> geoms(Nlev);
    //nGrow is how many cells we need outside of the box
    int nGrow = 0;
    //Loop over levels
    for (int lev=0; lev<Nlev; ++lev)
    {
      //Get the array of boxes for this level
      const BoxArray ba = amrData.boxArray(lev);
      //Distribution mapping i.e. how are boxes distributed across processors
      const DistributionMapping dm(ba);
      //set up a multifab for the outdata with the same boxarray, distribution mapping, components and ngrow
      outdata[lev] = MultiFab(ba,dm,nCompOut,nGrow);
      //we'll load each level as we go, so only need to fill one multifab at a time
      MultiFab indata(ba,dm,nCompIn,nGrow);
      int coord = 0; //indicates cartesian
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps); //magic IO call
      Print() << "Data has been read for level " << lev << std::endl;
      //Get array of arrays on this level - first array corresponds to the box, after that is ijk and component       
      auto const& inarrays  = indata.const_arrays();
      auto const& outarrays = outdata[lev].arrays();
      
      //If you can abstract to higher than the MPI process, you can do things here 
      ParallelFor(indata,[=] AMREX_GPU_DEVICE (int bx, int i, int j, int k) {
			   for (int n = 0; n < numProgressVars; n++) {
			     outarrays[bx](i,j,k,n) = (inarrays[bx](i,j,k,n)-d_c_u[n])/(d_c_b[n]-d_c_u[n]);
			   }
			 });
      
      //If you can't, use a loop like this.
      /*
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  auto const& out_a = outdata[lev].array(mfi);
	  auto const& in_a = indata.array(mfi);
	  ParallelFor(bx, [=]
	  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	  {
	  // Do something funky
	  for (int n = 0; n < nComp; n++) {
	  out_a(i,j,k,n) = in_a(i,j,k,n); 
	  }
	  });
	}
      */

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_prog");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev,0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2,2,2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata),outNames, geoms, 0.0, isteps, refRatios);
  }
  amrex::Finalize();
  return 0;
}
