#ifndef PARTICLEBASIS_H_
#define PARTICLEBASIS_H_

#include<constants.H>
#include <AMReX.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real hat1d(amrex::Real xi,amrex::Real xp,amrex::Real hatsize)
{
    amrex::Real zeta=xi-xp;
    amrex::Real halfhatsize=hatsize*half;
    amrex::Real funcval=zero;
    amrex::Real maxval=one/halfhatsize;
    amrex::Real slope=four/(hatsize*hatsize);

    if(amrex::Math::abs(zeta) < halfhatsize)
    {
        if(zeta>0.0)
        {
           funcval = -slope*zeta + maxval;       
        }
        else
        {
           funcval =  slope*zeta + maxval;
        }
    }
    return(funcval);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real hat3d(amrex::Real xi[AMREX_SPACEDIM], amrex::Real xp[AMREX_SPACEDIM],
                amrex::Real hatsize[AMREX_SPACEDIM])
{
    amrex::Real funcval=one;
    funcval *= hat1d(xi[XDIR],xp[XDIR],hatsize[XDIR]);
    funcval *= hat1d(xi[YDIR],xp[YDIR],hatsize[YDIR]);
    funcval *= hat1d(xi[ZDIR],xp[ZDIR],hatsize[ZDIR]);

    return(funcval);
}

#endif
