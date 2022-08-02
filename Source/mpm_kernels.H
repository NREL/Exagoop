#ifndef MPM_KERNELS_H_
#define MPM_KERNELS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <constants.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int  applybc(Real relvel_in[AMREX_SPACEDIM],
        Real relvel_out[AMREX_SPACEDIM],
        Real wall_mu, Real normaldir[AMREX_SPACEDIM], int bc)
{
    
    //look at section 8
    //Stomakhin et al., "A material point method for snow simulation",
    //ACS Trans. Graphics, 32, 4, 102, 2013

    relvel_out[XDIR]=relvel_in[XDIR];
    relvel_out[YDIR]=relvel_in[YDIR];
    relvel_out[ZDIR]=relvel_in[ZDIR];
    int modify_pos=false;

    if(bc==BC_SLIPWALL || bc==BC_NOSLIPWALL || bc==BC_PARTIALSLIPWALL)
    {
        Real velt[AMREX_SPACEDIM],velt_hat[AMREX_SPACEDIM];
        
        //normal is into the domain
        Real veln=relvel_in[XDIR]*normaldir[XDIR]+relvel_in[YDIR]*normaldir[YDIR]
            +relvel_in[ZDIR]*normaldir[ZDIR];

        velt[XDIR]=relvel_in[XDIR]-veln*normaldir[XDIR];
        velt[YDIR]=relvel_in[YDIR]-veln*normaldir[YDIR];
        velt[ZDIR]=relvel_in[ZDIR]-veln*normaldir[ZDIR];

        Real veltmag=std::sqrt(velt[XDIR]*velt[XDIR]
        +velt[YDIR]*velt[YDIR]+velt[ZDIR]*velt[ZDIR]);

        velt_hat[XDIR]=velt[XDIR]/(veltmag+TINYVAL);
        velt_hat[YDIR]=velt[YDIR]/(veltmag+TINYVAL);
        velt_hat[ZDIR]=velt[ZDIR]/(veltmag+TINYVAL);

        if(veln<zero)
        {
            modify_pos=true;

            if(bc==BC_SLIPWALL)
            {
                relvel_out[XDIR]=velt[XDIR];
                relvel_out[YDIR]=velt[YDIR];
                relvel_out[ZDIR]=velt[ZDIR];
            }
            else if(bc==BC_NOSLIPWALL)
            {
                relvel_out[XDIR]=zero;
                relvel_out[YDIR]=zero;
                relvel_out[ZDIR]=zero;
            }
            else if(bc==BC_PARTIALSLIPWALL)
            {
                if(veltmag <= -wall_mu*veln)
                {
                    relvel_out[XDIR]=zero;
                    relvel_out[YDIR]=zero;
                    relvel_out[ZDIR]=zero;
                }
                else
                {
                    relvel_out[XDIR]=(veltmag+wall_mu*veln)*velt_hat[XDIR];
                    relvel_out[YDIR]=(veltmag+wall_mu*veln)*velt_hat[YDIR];
                    relvel_out[ZDIR]=(veltmag+wall_mu*veln)*velt_hat[ZDIR];
                }
            }
            else
            {
                amrex::Abort("Unknown boundary condition\n");
            }
        }
    }
    return(modify_pos);
}

#endif
