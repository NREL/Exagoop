#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <constants.H>

using namespace amrex;

void write_plot_file(std::string fname, MultiFab &nodaldata, Vector<std::string> fieldnames, 
                           Geometry geom, BoxArray ba, DistributionMapping dm,Real time)
{
  MultiFab plotmf(ba,dm,nodaldata.nComp(),0);
  amrex::average_node_to_cellcenter(plotmf, 0, nodaldata, 0, nodaldata.nComp());
  WriteSingleLevelPlotfile(fname, plotmf, fieldnames, geom, time, 0);
}

void nodal_update(MultiFab &nodaldata,const amrex::Real& dt)
{
  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
  {
        const Box& bx=mfi.validbox();
        Box nodalbox = convert(bx, {1, 1, 1});
        
        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        amrex::ParallelFor(nodalbox,[=]
                AMREX_GPU_DEVICE (int i,int j,int k) noexcept
        {
           if(nodal_data_arr(i,j,k,MASS_INDEX) > zero)
           {
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    nodal_data_arr(i,j,k,VELX_INDEX+d) 
                    += nodal_data_arr(i,j,k,FRCX_INDEX+d)/nodal_data_arr(i,j,k,MASS_INDEX)*dt;
                }
           }  
        
        });
  }
}

void nodal_bcs(const amrex::Geometry geom,
               MultiFab &nodaldata,const amrex::Real& dt)
{
    //FIXME: This is a hard dirichlet for now
    
    const int* domloarr = geom.Domain().loVect();
    const int* domhiarr = geom.Domain().hiVect();

    GpuArray<int,AMREX_SPACEDIM> domlo={domloarr[0],domloarr[1],domloarr[2]};
    GpuArray<int,AMREX_SPACEDIM> domhi={domhiarr[0],domhiarr[1],domhiarr[2]};

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        const Box& bx=mfi.validbox();
        Box nodalbox = convert(bx, {1, 1, 1});

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        amrex::ParallelFor(nodalbox,[=]
        AMREX_GPU_DEVICE (int i,int j,int k) noexcept
        {
            IntVect nodeid(i,j,k);

            if(   nodeid[XDIR]==domlo[XDIR] || nodeid[XDIR]==(domhi[XDIR]+1)
              ||  nodeid[YDIR]==domlo[YDIR] || nodeid[YDIR]==(domhi[YDIR]+1)
              ||  nodeid[ZDIR]==domlo[ZDIR] || nodeid[ZDIR]==(domhi[ZDIR]+1) )
            {
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    nodal_data_arr(i,j,k,VELX_INDEX+d)=zero;
                }
            }  

        });
    }
}
