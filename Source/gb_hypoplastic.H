#ifndef GB_HYPOPLASTIC_H_
#define GB_HYPOPLASTIC_H_

#include<constants.H>
#include <AMReX.H>
#include<mpm_particle_container.H>

/**
 * @brief Yudong: implementation of Gudehus-Bauer hypoplastic model 
 * 
 * @param epsdot 
 * @param sigma 
 * @param e : void ratio for material properties
 * @param dt : time step for updating stress tensor
 * @return AMREX_GPU_DEVICE 
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void GB_hypoplastic(amrex::Real eps_dot[NCOMP_TENSOR],
        amrex::Real omega_dot[NCOMP_TENSOR],
        amrex::Real sigma[NCOMP_TENSOR],
        amrex::Real dt,
        MPMParticleContainer::ParticleType &p,
        MPMspecs specs
        )
{
    
    // normalized stress tensor
    amrex::Real sigma_hat[NCOMP_TENSOR]; 
    // normalized deviatoric stress tenstor
    amrex::Real sigma_hat_star[NCOMP_TENSOR];
    // Jaumann stress rate
    amrex::Real sigma_o[NCOMP_TENSOR];
    // stress rate tensor
    amrex::Real sigma_dot[NCOMP_TENSOR];
    // identity matrix upper
    amrex::Real Identity[NCOMP_TENSOR] = {1.0,0,0,1.0,0,1.0};
    // void ratio
    amrex::Real e = p.rdata(realData::void_ratio);

    //get tensors for computation of stress rate
    amrex::Real tensor_Jaumann[AMREX_SPACEDIM*AMREX_SPACEDIM];
    amrex::Real tensor_spin[AMREX_SPACEDIM*AMREX_SPACEDIM];
    amrex::Real tensor_stress[AMREX_SPACEDIM*AMREX_SPACEDIM];

    // intermediate products
    amrex::Real tensor_omega_dot_sigma[AMREX_SPACEDIM*AMREX_SPACEDIM];
    amrex::Real tensor_sigma_dot_omega[AMREX_SPACEDIM*AMREX_SPACEDIM];

    // Linear contribution
    amrex::Real L[NCOMP_TENSOR];
    // Nonlinear contribution
    amrex::Real N[NCOMP_TENSOR];

    // GH hypoplastic model parameters
    amrex::Real phi_c = specs.GBparameters.phi_c; // critical state internal friction angle
    amrex::Real h_s = specs.GBparameters.h_s; // granular hardness
    amrex::Real n = specs.GBparameters.n; // pressure sensitive exponent : G paper says range (0.3, 0.5)
    amrex::Real e_d0 = specs.GBparameters.e_d0; // material constant for minimum void ratio
    amrex::Real e_c0 = specs.GBparameters.e_c0; // material constant for critical void ratio
    amrex::Real e_i0 = specs.GBparameters.e_i0; // material constant for maximum void ratio
    amrex::Real alpha = specs.GBparameters.alpha; //material constant (0.1, 0.3)  
    amrex::Real beta = specs.GBparameters.beta; //material constant (1.0, 1.1)
    amrex::Real stress_cutoff = specs.GBparameters.stress_cutoff;

    amrex::Real trace_sigma = sigma[XX] + sigma[YY] + sigma[ZZ];

    if(trace_sigma>0){
        trace_sigma = 0;
        if(e>e_i0){
            e = e_i0;
        }
        else if (e<e_d0){
            e = e_d0;
        }
        p.rdata(realData::void_ratio) = e;
    }

    if(trace_sigma != 0.0){
        // compute some varaibles
        amrex::Real trace_eps_dot=eps_dot[XX]+eps_dot[YY]+eps_dot[ZZ];
        amrex::Real modulus_eps_dot = sqrt(pow(trace_eps_dot,2));

        for(int i = 0; i < NCOMP_TENSOR; i++) {
            sigma_hat[i] = sigma[i]/trace_sigma;
            sigma_hat_star[i] = sigma_hat[i] - (1.0/3.0)*Identity[i];
        }
        amrex::Real trace_sigma_hat_star = sigma_hat_star[XX]+sigma_hat_star[YY]+sigma_hat_star[ZZ];
        //amrex::Real modulus_sigma = sqrt(pow(trace_sigma,2));
            

        amrex::Real modulus_sigma_hat_star = sqrt(pow(trace_sigma_hat_star,2));

        amrex::Real trace_sigma_hat_star_square = pow(sigma_hat_star[XX],2)+pow(sigma_hat_star[YY],2)+pow(sigma_hat_star[ZZ],2)
                                                    +2*pow(sigma_hat_star[XY],2)+2*pow(sigma_hat_star[YZ],2)+2*pow(sigma_hat_star[XZ],2);
        amrex::Real det_sigma_hat_star = sigma_hat_star[XX]*sigma_hat_star[YY]*sigma_hat_star[ZZ]
                                        +2*sigma_hat_star[XY]*sigma_hat_star[YZ]*sigma_hat_star[XZ]
                                        -sigma_hat_star[XY]*sigma_hat_star[XY]*sigma_hat_star[ZZ]
                                        -sigma_hat_star[YZ]*sigma_hat_star[YZ]*sigma_hat_star[XX]
                                        -sigma_hat_star[XZ]*sigma_hat_star[XZ]*sigma_hat_star[YY];
        amrex::Real trace_sigma_hat_star_cube = 3*det_sigma_hat_star;

        amrex::Real trace_sigma_hat_dot_eps_dot = sigma_hat[XX]*eps_dot[XX] + sigma_hat[YY]*eps_dot[YY] +sigma_hat[ZZ]*eps_dot[ZZ]
                                            + 2*sigma_hat[XY]*eps_dot[XY] + 2*sigma_hat[YZ]*eps_dot[YZ] + 2*sigma_hat[XZ]*eps_dot[XZ];
        
        // compute c1 c2
        amrex::Real c_1 = sqrt(3.0/8.0)*((3-sin(phi_c))/sin(phi_c));
        amrex::Real c_2 = (3.0/8.0)*((3+sin(phi_c))/sin(phi_c));

        // update stiffness factor
        amrex::Real h_i = 1/pow(c_1,2) + (1.0/3.0)*pow((e_i0+e_d0)/(e_c0-e_d0),alpha)*(1/(sqrt(3)*c_1));

        // update void ratio related variables
        amrex::Real e_i = e_i0*exp(-pow(-trace_sigma/h_s,n)); // max void ratio
        amrex::Real e_d = e_d0*exp(-pow(-trace_sigma/h_s,n)); // min void ratio
        amrex::Real e_c = e_c0*exp(-pow(-trace_sigma/h_s,n)); // critical void ratio
        p.rdata(realData::min_void_ratio) = e_d;
        p.rdata(realData::max_void_ratio) = e_i;
        p.rdata(realData::critical_void_ratio) = e_c;

        if(e>e_i){
            e = e_i;
        }
        else if (e<e_d){
            e = e_d;
        }

        // Pycnotropy factors
        // amrex::Real f_d = pow((e-e_d)/(e_c-e_d),alpha);
        amrex::Real f_d = (e-e_d)/(e_c-e_d);
        if(f_d>0){
            f_d = pow(f_d,alpha);
        } else {
            f_d = 0;
        }
        amrex::Real f_e = pow(e_c/e, beta);
        amrex::Real f_b = pow(e_i0/e_c0, beta)*(h_s/(n*h_i))*((1+e_i)/e_i)*pow(-trace_sigma/h_s,1-n);
        amrex::Real f_s = f_b*f_e;

        // compute lode angle
        // amrex::Real cos3theta = -sqrt(6.0)*trace_sigma_hat_star_cube/pow(trace_sigma_hat_star_square,1.5);
        amrex::Real cos3theta;
        if(trace_sigma_hat_star_square<=1e-10){
            cos3theta = 1.0;
        } else {
            cos3theta = -sqrt(6.0)*trace_sigma_hat_star_cube/pow(trace_sigma_hat_star_square,1.5);
            if(cos3theta>1.0) cos3theta = 1.0;
            if(cos3theta<-1.0) cos3theta = -1.0;
        }

        amrex::Real a_1 = 1.0/(c_1+c_2*modulus_sigma_hat_star*(1+cos3theta));
        
        // update stress based on stress rate
        for(int i=0; i<NCOMP_TENSOR; i++){
            N[i] = a_1*(sigma_hat[i]+sigma_hat_star[i]);
            // linear contribution
            L[i] = pow(a_1,2)*eps_dot[i]+sigma_hat[i]*trace_sigma_hat_dot_eps_dot;
            sigma_o[i] = f_s*(L[i]+f_d*N[i]*modulus_eps_dot);
        }
        
        get_antisymmetric_tensor(omega_dot,tensor_spin);
        get_tensor_from_array(sigma, tensor_stress);
        get_tensor_from_array(sigma_o, tensor_Jaumann);
        
        tensor_tensor_dot_product(tensor_spin,tensor_stress,tensor_omega_dot_sigma);
        tensor_tensor_dot_product(tensor_stress,tensor_spin,tensor_sigma_dot_omega);
        
        int index = 0;
        for(int i=0; i<AMREX_SPACEDIM; i++){
            for(int j=i; j<AMREX_SPACEDIM; j++){
                sigma_dot[index] = tensor_Jaumann[i*AMREX_SPACEDIM+j]
                                 + tensor_omega_dot_sigma[i*AMREX_SPACEDIM+j]
                                 - tensor_sigma_dot_omega[i*AMREX_SPACEDIM+j];
                sigma[index] = sigma[index]+sigma_dot[index]*dt;
                index++;
            }   
        }
        //limit stress to be compressive
        if(sigma[XX]>0){
            sigma[XX] = -10.0;
            sigma[XY] = 0.0;
            sigma[XZ] = 0.0;
        }
        if(sigma[YY]>0){
            sigma[YY] = -10.0;
            sigma[YZ] = 0.0;
        }
        if(sigma[ZZ]>0){
            sigma[ZZ] = -10.0;
        }

        // update void ratio here
        amrex::Real old_void_ratio = e;
        e = old_void_ratio + (p.rdata(realData::strainrate+XX)
                                    +p.rdata(realData::strainrate+YY)
                                    +p.rdata(realData::strainrate+ZZ))
                                    * (1 +old_void_ratio)* dt;

        trace_sigma = sigma[XX] + sigma[YY] + sigma[ZZ];
        e_i = e_i0*exp(-pow(-trace_sigma/h_s,n)); // max void ratio
        e_d = e_d0*exp(-pow(-trace_sigma/h_s,n)); // min void ratio
        e_c = e_c0*exp(-pow(-trace_sigma/h_s,n)); // critical void ratio
        p.rdata(realData::min_void_ratio) = e_d;
        p.rdata(realData::max_void_ratio) = e_i;
        p.rdata(realData::critical_void_ratio) = e_c;
        if(e>e_i){
            e = e_i;
        }
        else if (e<e_d){
            e = e_d;
        }
        p.rdata(realData::void_ratio) = e;
        p.rdata(realData::volume)	= (1+e)/(1+old_void_ratio)*p.rdata(realData::vol_init); // use order_of_scheme=1 to update vol_init
        p.rdata(realData::density)	= p.rdata(realData::mass)/p.rdata(realData::volume);

        //amrex::Print()<<"strain rate XX = " << p.rdata(realData::strainrate+XX);
        //e += (p.rdata(realData::strainrate+XX)+p.rdata(realData::strainrate+YY) +p.rdata(realData::strainrate+ZZ)) * (1 +e)* dt;
        //sigma[XX] = (sigma[XX]<=stress_cutoff)?sigma[XX]:stress_cutoff;
        //sigma[YY] = (sigma[YY]<=stress_cutoff)?sigma[YY]:stress_cutoff;
        //sigma[ZZ] = (sigma[ZZ]<=stress_cutoff)?sigma[ZZ]:stress_cutoff;
    }
    else{

        if(e>e_i0){
            e = e_i0;
        }
        else if (e<e_d0){
            e = e_d0;
        }
        p.rdata(realData::void_ratio) = e;
    }
}


#endif
