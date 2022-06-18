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

void backup_current_velocity(MultiFab &nodaldata)
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
        	   	nodal_data_arr(i,j,k,MASS_OLD_INDEX) = nodal_data_arr(i,j,k,MASS_INDEX);
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                	//Storing V(t) values at these memory locations. This is to be followed by a call to store_delta_velocity
                    nodal_data_arr(i,j,k,DELTA_VELX_INDEX+d) = nodal_data_arr(i,j,k,VELX_INDEX+d);
                }
           }

        });
  }
}



void store_delta_velocity(MultiFab &nodaldata)
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
                    nodal_data_arr(i,j,k,DELTA_VELX_INDEX+d) = nodal_data_arr(i,j,k,VELX_INDEX+d)-nodal_data_arr(i,j,k,DELTA_VELX_INDEX+d);
                }
           }

        });
  }
}

void nodal_update(MultiFab &nodaldata,const amrex::Real& dt, const amrex::Real& mass_tolerance)
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
                	if(nodal_data_arr(i,j,k,MASS_INDEX)>=mass_tolerance)
                	{
                		nodal_data_arr(i,j,k,VELX_INDEX+d) += nodal_data_arr(i,j,k,FRCX_INDEX+d)/nodal_data_arr(i,j,k,MASS_INDEX)*dt;
                	}
                	else
                	{
                		nodal_data_arr(i,j,k,VELX_INDEX+d) = 0.0;
                	}

                }
           }
        
        });
  }
}

void initialise_shape_function_indices(iMultiFab &shapefunctionindex,const amrex::Geometry geom)
{
	const int* domloarr = geom.Domain().loVect();
	const int* domhiarr = geom.Domain().hiVect();

    int periodic[AMREX_SPACEDIM]={geom.isPeriodic(XDIR),
	                                  geom.isPeriodic(YDIR),
	                                  geom.isPeriodic(ZDIR)};

    GpuArray<int,AMREX_SPACEDIM> domlo={domloarr[0],domloarr[1],domloarr[2]};
    GpuArray<int,AMREX_SPACEDIM> domhi={domhiarr[0],domhiarr[1],domhiarr[2]};

    const auto domain = geom.Domain();
    const int* lo = domain.loVect ();
    const int* hi = domain.hiVect ();

    for (MFIter mfi(shapefunctionindex); mfi.isValid(); ++mfi)
    {
        const Box& bx=mfi.validbox();
        Box nodalbox = convert(bx, {1, 1, 1});

        Array4<int> shapefunctionindex_arr=shapefunctionindex.array(mfi);

	    amrex::ParallelFor(nodalbox,[=]
	    AMREX_GPU_DEVICE (int i,int j,int k) noexcept
	    {
	    	if(i==lo[0])
    		{
	    		shapefunctionindex_arr(i,j,k,0)=0;
    		}
	    	else if(i==lo[0]+1)
	    	{
	    		shapefunctionindex_arr(i,j,k,0)=1;
	    	}
	    	else if(i==hi[0])
	    	{
	    		shapefunctionindex_arr(i,j,k,0)=4;
	    	}
	    	else if(i==hi[0]-1)
	    	{
	    		shapefunctionindex_arr(i,j,k,0)=3;
	    	}
	    	else
	    	{
	    		shapefunctionindex_arr(i,j,k,0)=2;
	    	}

	    	if(j==lo[1])
	    	{
	    		shapefunctionindex_arr(i,j,k,1)=0;
	    	}
	    	else if(j==lo[1]+1)
	    	{
	    		shapefunctionindex_arr(i,j,k,1)=1;
	    	}
	    	else if(j==hi[1])
	    	{
	    		shapefunctionindex_arr(i,j,k,1)=4;
	    	}
	    	else if(j==hi[1]-1)
	    	{
	    		shapefunctionindex_arr(i,j,k,1)=3;
	    	}
	    	else
	    	{
	    		shapefunctionindex_arr(i,j,k,1)=2;
	    	}

	    	if(k==lo[2])
	    	{
	    		shapefunctionindex_arr(i,j,k,2)=0;
	    	}
	    	else if(k==lo[2]+1)
	    	{
	    		shapefunctionindex_arr(i,j,k,2)=1;
	    	}
	    	else if(k==hi[2])
	    	{
	    		shapefunctionindex_arr(i,j,k,2)=4;
	    	}
	    	else if(k==hi[2]-1)
	    	{
	    		shapefunctionindex_arr(i,j,k,2)=3;
	    	}
	    	else
	    	{
	    		shapefunctionindex_arr(i,j,k,2)=2;
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
    
    int periodic[AMREX_SPACEDIM]={geom.isPeriodic(XDIR),
                                  geom.isPeriodic(YDIR),
                                  geom.isPeriodic(ZDIR)};

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
            bool impose_wall_x=false;
            bool impose_wall_y=false;
            bool impose_wall_z=false;

            impose_wall_x=!periodic[XDIR] && 
                        (nodeid[XDIR]==domlo[XDIR] || nodeid[XDIR]==(domhi[XDIR]+1))?true:false;
            
            impose_wall_y=!periodic[YDIR] && 
                        (nodeid[YDIR]==domlo[YDIR] || nodeid[YDIR]==(domhi[YDIR]+1))?true:false;

            impose_wall_z=!periodic[ZDIR] && 
                        (nodeid[ZDIR]==domlo[ZDIR] || nodeid[ZDIR]==(domhi[ZDIR]+1))?true:false;


            if(impose_wall_x)
            {
            	nodal_data_arr(i,j,k,VELX_INDEX+0)=zero;
            	nodal_data_arr(i,j,k,VELX_INDEX+1)=zero;
            	nodal_data_arr(i,j,k,VELX_INDEX+2)=zero;
            }

            if(impose_wall_y)
            {
            	//nodal_data_arr(i,j,k,VELX_INDEX+0)=zero;
            	nodal_data_arr(i,j,k,VELX_INDEX+1)=zero;
            	//nodal_data_arr(i,j,k,VELX_INDEX+2)=zero;
            }

            if(impose_wall_z)
            {
            	nodal_data_arr(i,j,k,VELX_INDEX+2)=zero;
            }

            /*
            if(impose_wall_x || impose_wall_y || impose_wall_z)
            {
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    nodal_data_arr(i,j,k,VELX_INDEX+d)=zero;
                }
            } */

        });
    }
}
