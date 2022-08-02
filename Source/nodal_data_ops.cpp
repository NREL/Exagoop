#include <nodal_data_ops.H>
#include <mpm_eb.H>
#include <mpm_kernels.H>

using namespace amrex;

void write_plot_file(std::string fname, MultiFab &nodaldata, Vector<std::string> fieldnames, 
                           Geometry geom, BoxArray ba, DistributionMapping dm,Real time)
{
  MultiFab plotmf(ba,dm,nodaldata.nComp(),0);
  average_node_to_cellcenter(plotmf, 0, nodaldata, 0, nodaldata.nComp());
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

void nodal_levelset_bcs(MultiFab &nodaldata,const Geometry geom,amrex::Real &dt,int slip)
{
  //need something more sophisticated
  //but lets get it working!
  //
    int lsref=mpm_ebtools::ls_refinement;
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dx = geom.CellSizeArray();

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        const Box& bx=mfi.validbox();
        Box nodalbox = convert(bx, {1, 1, 1});

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);
        Array4<Real> lsarr=mpm_ebtools::lsphi->array(mfi);

        amrex::ParallelFor(nodalbox,[=]
        AMREX_GPU_DEVICE (int i,int j,int k) noexcept
        {
            IntVect nodeid(i,j,k);
            IntVect refined_nodeid(i*lsref,j*lsref,k*lsref);

            if(!slip)
            {

                if(lsarr(refined_nodeid) < TINYVAL)
                {
                    nodal_data_arr(nodeid,VELX_INDEX+XDIR)=zero;
                    nodal_data_arr(nodeid,VELX_INDEX+YDIR)=zero;
                    nodal_data_arr(nodeid,VELX_INDEX+ZDIR)=zero;
                }
            }
            else
            {
                if(lsarr(refined_nodeid) < TINYVAL)
                {
                    amrex::Real eps=0.00001;

                    amrex::Real xp[AMREX_SPACEDIM]={plo[XDIR]+(i+eps)*dx[XDIR],
                        plo[YDIR]+(j+eps)*dx[YDIR],plo[ZDIR]+(k+eps)*dx[ZDIR]};

                    amrex::Real norm[AMREX_SPACEDIM]={1.0,0.0,0.0};
                    amrex::Real dist=get_levelset_value(lsarr,plo,dx,xp,lsref); 
                    dist=amrex::Math::abs(dist);

                    get_levelset_grad(lsarr,plo,dx,xp,lsref,norm);
                    amrex::Real gradmag=std::sqrt(norm[XDIR]*norm[XDIR]
                                                + norm[YDIR]*norm[YDIR]
                                                + norm[ZDIR]*norm[ZDIR]);

                    if(gradmag>eps)
                    {
                        norm[XDIR]=norm[XDIR]/gradmag;
                        norm[YDIR]=norm[YDIR]/gradmag;
                        norm[ZDIR]=norm[ZDIR]/gradmag;

                        amrex::Real reflect_point[AMREX_SPACEDIM]=
                        { xp[XDIR]+two*dist*norm[XDIR], 
                          xp[YDIR]+two*dist*norm[YDIR], 
                          xp[ZDIR]+two*dist*norm[ZDIR] };

                        amrex::Real sum_dist2=0.0;
                        amrex::Real velx_dist2=0.0;
                        amrex::Real vely_dist2=0.0;
                        amrex::Real velz_dist2=0.0;

                        for(int kk=-1;kk<1;kk++)
                        {
                            for(int jj=-1;jj<1;jj++)
                            {
                                for(int ii=-1;ii<1;ii++)
                                {
                                    if(ii!=0 && jj!=0 && kk!=0)
                                    {

                                        if(lsarr((i+ii)*lsref,(j+jj)*lsref,(k+kk)*lsref) > TINYVAL)
                                        {
                                            amrex::Real xp[AMREX_SPACEDIM]=
                                            { plo[XDIR]+(i+ii)*dx[XDIR], 
                                              plo[YDIR]+(j+jj)*dx[YDIR], 
                                              plo[ZDIR]+(k+kk)*dx[ZDIR] };

                                            amrex::Real dist2 = std::pow((xp[XDIR]-reflect_point[XDIR]),2.0) +
                                                                std::pow((xp[YDIR]-reflect_point[YDIR]),2.0) +
                                                                std::pow((xp[ZDIR]-reflect_point[ZDIR]),2.0);

                                            sum_dist2 += 1.0/dist2;

                                            velx_dist2 += nodal_data_arr(i+ii,j+jj,k+kk,VELX_INDEX)/dist2; 
                                            vely_dist2 += nodal_data_arr(i+ii,j+jj,k+kk,VELY_INDEX)/dist2; 
                                            velz_dist2 += nodal_data_arr(i+ii,j+jj,k+kk,VELZ_INDEX)/dist2; 

                                        }
                                    }
                                }
                            }
                        }

                        if(sum_dist2 > 0.0)
                        {
                            amrex::Real vel_refl_pnt[AMREX_SPACEDIM]=
                            {velx_dist2/sum_dist2, vely_dist2/sum_dist2, velz_dist2/sum_dist2};

                            amrex::Real vel_n = vel_refl_pnt[XDIR]*norm[XDIR]
                                              + vel_refl_pnt[YDIR]*norm[YDIR]
                                              + vel_refl_pnt[ZDIR]*norm[ZDIR];

                            amrex::Real velt[AMREX_SPACEDIM]=
                            {vel_refl_pnt[XDIR]-vel_n*norm[XDIR],
                             vel_refl_pnt[YDIR]-vel_n*norm[YDIR],
                             vel_refl_pnt[ZDIR]-vel_n*norm[ZDIR]};


                            nodal_data_arr(i,j,k,VELX_INDEX)=-vel_n*norm[XDIR]+velt[XDIR];
                            nodal_data_arr(i,j,k,VELY_INDEX)=-vel_n*norm[YDIR]+velt[YDIR];
                            nodal_data_arr(i,j,k,VELZ_INDEX)=-vel_n*norm[ZDIR]+velt[ZDIR];
                        }
                        else
                        {
                            nodal_data_arr(i,j,k,VELX_INDEX)=0.0;
                            nodal_data_arr(i,j,k,VELY_INDEX)=0.0;
                            nodal_data_arr(i,j,k,VELZ_INDEX)=0.0;
                        }

                    }


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
                    nodal_data_arr(i,j,k,DELTA_VELX_INDEX+d) 
                    = nodal_data_arr(i,j,k,VELX_INDEX+d)-nodal_data_arr(i,j,k,DELTA_VELX_INDEX+d);
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
            if(nodal_data_arr(i,j,k,MASS_INDEX) >=mass_tolerance)
            {
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                        nodal_data_arr(i,j,k,VELX_INDEX+d) += 
                        nodal_data_arr(i,j,k,FRCX_INDEX+d)/nodal_data_arr(i,j,k,MASS_INDEX)*dt;
                }
            }
            else
            {
                nodal_data_arr(i,j,k,VELX_INDEX) = 0.0;
                nodal_data_arr(i,j,k,VELY_INDEX) = 0.0;
                nodal_data_arr(i,j,k,VELZ_INDEX) = 0.0;
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
        MultiFab &nodaldata,int bclo[AMREX_SPACEDIM],
        int bchi[AMREX_SPACEDIM],Real wall_mu_lo[AMREX_SPACEDIM],
        Real wall_mu_hi[AMREX_SPACEDIM],const amrex::Real& dt)
{

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
            Real relvel_in[AMREX_SPACEDIM],relvel_out[AMREX_SPACEDIM];

            for(int d=0;d<AMREX_SPACEDIM;d++)
            {
                relvel_in[d]=nodal_data_arr(nodeid,VELX_INDEX+d);
                relvel_out[d]=relvel_in[d];
            }

            if(nodeid[XDIR]==domlo[XDIR])
            {
                Real normaldir[AMREX_SPACEDIM]={1.0,0.0,0.0};
                int tmp=applybc(relvel_in,relvel_out,wall_mu_lo[XDIR],
                        normaldir,bclo[XDIR]);
            }
            if(nodeid[XDIR]==(domhi[XDIR]+1))
            {
                Real normaldir[AMREX_SPACEDIM]={-1.0,0.0,0.0};
                int tmp=applybc(relvel_in,relvel_out,wall_mu_hi[XDIR],
                        normaldir,bchi[XDIR]);
            }
            if(nodeid[YDIR]==domlo[YDIR])
            {
                Real normaldir[AMREX_SPACEDIM]={0.0,1.0,0.0};
                int tmp=applybc(relvel_in,relvel_out,wall_mu_lo[YDIR],
                        normaldir,bclo[YDIR]);
            }
            if(nodeid[YDIR]==(domhi[YDIR]+1))
            {
                Real normaldir[AMREX_SPACEDIM]={0.0,-1.0,0.0};
                int tmp=applybc(relvel_in,relvel_out,wall_mu_hi[YDIR],
                        normaldir,bchi[YDIR]);
            }
            if(nodeid[ZDIR]==domlo[ZDIR])
            {
                Real normaldir[AMREX_SPACEDIM]={0.0,0.0,1.0};
                int tmp=applybc(relvel_in,relvel_out,wall_mu_lo[ZDIR],
                        normaldir,bclo[ZDIR]);
            }
            if(nodeid[ZDIR]==(domhi[ZDIR]+1))
            {
                Real normaldir[AMREX_SPACEDIM]={0.0,0.0,-1.0};
                int tmp=applybc(relvel_in,relvel_out,wall_mu_hi[ZDIR],
                        normaldir,bchi[ZDIR]);
            }
            
            for(int d=0;d<AMREX_SPACEDIM;d++)
            {
                nodal_data_arr(nodeid,VELX_INDEX+d)=relvel_out[d];
            }

        });
    }
}
