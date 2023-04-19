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
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void miuI(amrex::Real eps_dot[NCOMP_TENSOR],
        amrex::Real sigma[NCOMP_TENSOR],
        MPMParticleContainer::ParticleType &p,
        MPMspecs specs
        )
{
    //input
    //amrex::Print()<<"now in miuI model";
    amrex::Real d = specs.miu_I_Parameters.d;
    amrex::Real rou_s = specs.miu_I_Parameters.rou_s;
    amrex::Real miu_s = specs.miu_I_Parameters.miu_s;
    amrex::Real miu_2 = specs.miu_I_Parameters.miu_2;
    amrex::Real I_0 = specs.miu_I_Parameters.I_0;

    amrex::Real tau[NCOMP_TENSOR] = {0.0,0,0,0.0,0,0.0}; 
    amrex::Real delta[NCOMP_TENSOR] = {1.0,0,0,1.0,0,1.0};

    amrex::Real one_by_three=1.0/3.0;
    amrex::Real trace_sigma = sigma[XX]+sigma[YY]+sigma[ZZ];
    amrex::Real P = - one_by_three * trace_sigma; //should be isostatic pressure?
    for(int i=0; i<NCOMP_TENSOR; i++){
        if(eps_dot[i]<TINYVAL){
            eps_dot[i] = 0.0;
        }
    }
    if(P<TINYVAL && P>0) {
        P = TINYVAL; 
    }
    else if (P>-TINYVAL && P<0){
        P = -TINYVAL;
    }
    else if (P == 0.0) {
        P = TINYVAL;
    }
    amrex::Real eps_dot_2nd_invariant = eps_dot[XX]*eps_dot[YY] + eps_dot[YY]*eps_dot[ZZ] + eps_dot[XX]*eps_dot[ZZ]
                                        - eps_dot[XY]*eps_dot[XY] - eps_dot[YZ]*eps_dot[YZ] - eps_dot[XZ]*eps_dot[XZ];
    if(eps_dot_2nd_invariant<TINYVAL && eps_dot_2nd_invariant>0) {
        eps_dot_2nd_invariant = TINYVAL; 
    }
    else if (eps_dot_2nd_invariant>-TINYVAL && eps_dot_2nd_invariant<0){
        eps_dot_2nd_invariant = TINYVAL;
    }
    else if (eps_dot_2nd_invariant == 0.0) {
        eps_dot_2nd_invariant = TINYVAL;
    }
    amrex::Real I = 0.0;
    // amrex::Print()<<"Pressure: " << P << "\n";
    // amrex::Print()<<"epsdot 2nd invariant: "<< eps_dot_2nd_invariant <<"\n";
    // amrex::Print()<<"rou_s: "<< rou_s <<"\n";
    // amrex::Print()<<"d: "<< d <<"\n";
    // amrex::Print()<<"denominator: sqrt P/rou_s : "<< pow(P/rou_s,0.5) <<"\n";
    I = eps_dot_2nd_invariant*d/pow(P/rou_s,0.5);
    p.rdata(realData::inertial_number) = I;
    if(I<TINYVAL && I>0) {
        I = TINYVAL; 
    }
    else if (I>-TINYVAL && I<0){
        I = TINYVAL;
    }
    else if (I == 0.0) {
        I = TINYVAL;
    }
    // amrex::Print()<<"Inertial number: "<< I << "\n";
    
    amrex::Real miu_I = miu_s + (miu_2 - miu_s)/(I_0/I + 1);
    p.rdata(realData::miu_I) = miu_I;
    // amrex::Print()<<"miuI: "<< miu_I << "\n";
    amrex::Real eta_eps_P = miu_I*P/eps_dot_2nd_invariant;
    // amrex::Print()<<"eta_eps_P: "<< eta_eps_P << "\n";
    for(int i=0; i<NCOMP_TENSOR; i++){
        tau[i] = eta_eps_P*eps_dot[i];
        sigma[i] = -P*delta[i] + tau[i];
        // amrex::Print()<<"tau["<<i<<"]: "<< tau[i] << "\n";
        // amrex::Print()<<"sigma["<<i<<"]: "<< sigma[i] << "\n";
    }
    for(int i=0; i<NCOMP_TENSOR; i++){
        p.rdata(realData::tau+i) = tau[i];
    }  
}
#endif
