#ifndef CONSTITUTIVE_MODELS_H_
#define CONSTITUTIVE_MODELS_H_

#include<constants.H>
#include <AMReX.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void linear_elastic(amrex::Real eps[AMREX_SPACEDIM*2],
        amrex::Real epsdot[AMREX_SPACEDIM*2],
        amrex::Real sigma[AMREX_SPACEDIM*2],
        amrex::Real E, amrex::Real v)
{
    amrex::Real trace_eps=eps[XX]+eps[YY]+eps[ZZ];
    amrex::Real const1=E/(one+v)/(one-two*v);
    amrex::Real const2=E/(one+v);
    
    sigma[XX]=const1*((1-v)*eps[XX] + v*(eps[YY]+eps[ZZ]));
    sigma[YY]=const1*((1-v)*eps[YY] + v*(eps[XX]+eps[ZZ]));
    sigma[ZZ]=const1*((1-v)*eps[ZZ] + v*(eps[XX]+eps[YY]));

    sigma[XY]=const2*eps[XY];
    sigma[XZ]=const2*eps[XZ];
    sigma[YZ]=const2*eps[YZ];
}

#endif
