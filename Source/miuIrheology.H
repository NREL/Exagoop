#ifndef MIUIRHEOLOGY_H_
#define MIUIRHEOLOGY_H_

#include<constants.H>
#include <AMReX.H>
#include<math.h>
#include<mpm_particle_container.H>

/**
 * @brief miuI rheology model implementation
 *  Yudong Li 01/05/2023
 * 
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void miuI(amrex::Real eps_dot[NCOMP_TENSOR],
        amrex::Real sigma[NCOMP_TENSOR],
        MPMParticleContainer::ParticleType &p,
        MPMspecs specs
        )
{
    //input
    amrex::Real d = specs.miu_I_Parameters.d;
    amrex::Real rou_s = specs.miu_I_Parameters.rou_s;
    amrex::Real miu_s = specs.miu_I_Parameters.miu_s;
    amrex::Real miu_2 = specs.miu_I_Parameters.miu_2;
    amrex::Real I_0 = specs.miu_I_Parameters.I_0;

    amrex::Real tau[NCOMP_TENSOR] = {0.0,0,0,0.0,0,0.0}; 
    amrex::Real delta[NCOMP_TENSOR] = {1.0,0,0,1.0,0,1.0};

    
    amrex::Real trace_sigma = sigma[XX]+sigma[YY]+sigma[ZZ];
    amrex::Real P = (1/3) * trace_sigma; //should be isostatic pressure?
    amrex::Real eps_dot_2nd_invariant = eps_dot[XX]*eps_dot[YY] + eps_dot[YY]*eps_dot[ZZ] + eps_dot[XX]*eps_dot[ZZ]
                                        - eps_dot[XY]*eps_dot[XY] - eps_dot[YZ]*eps_dot[YZ] - eps_dot[XZ]*eps_dot[XZ];
    amrex::Real I = eps_dot_2nd_invariant*d/pow(P/rou_s,0.5);
    p.rdata(realData::inertial_number) = I;
    if(I!=0){
        amrex::Real miu_I = miu_s + (miu_2 - miu_s)/(I_0/I + 1);
        p.rdata(realData::miu_I) = miu_I;
        amrex::Real eta_eps_P = miu_I*P/eps_dot_2nd_invariant;
        
        for(int i=0; i<NCOMP_TENSOR; i++){
            tau[i] = eta_eps_P*eps_dot[i];
            sigma[i] = -P*delta[i] + tau[i];
        }
    }
    for(int i=0; i<NCOMP_TENSOR; i++){
        p.rdata(realData::tau+i) = tau[i];
    }  
}
#endif
