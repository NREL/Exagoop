#pragma once
#include<math.h>
#include<iostream>
#include "constants.H"
#include "interpolants.H"

struct GBHypoParameters{
            //double phi_c = 0.5236; // critical state internal friction angle, unit: radian
            //double h_s = 5800000000.0; // granular hardness, unit: Pa
            //double n = 0.28; // pressure sensitive exponent : G paper says range (0.3, 0.5)
            //double e_d0 = 0.53; // material constant for minimum void ratio
            //double e_c0 = 0.84; // material constant for critical void ratio
            //double e_i0 = 1.0; // material constant for maximum void ratio
            //double alpha = 0.13; //material constant (0.1, 0.3)  
            //double beta = 1.0; //material constant (1.0, 1.1)

            // quartz sand properties
            double phi_c = 0.523599;
            double h_s = 190000000.0;
            double n = 0.4;
            double e_d0 = 0.51;
            double e_c0 = 0.82;
            double e_i0 = 1.2;
            double alpha = 0.15;
            double beta = 1.05;
            double stress_cutoff = -10.0; //stress cutoff to ensure compressive stress on material points
        };
/**
 * @brief Hypoplastic model
 * 
 * @param eps_dot  This is the strain rate tensor
 * @param omega_dot spin rate tensor (for a single element oedometer test, this is ignored)
 * @param sigma stress tensor
 * @param dt time step
 * @param e void ratio
 * @param GBparameters parameters for GB hypoplastic model 
 */
void GB_hypoplastic(double eps_dot[NCOMP_TENSOR],
        double omega_dot[NCOMP_TENSOR],
        double sigma[NCOMP_TENSOR],
        double dt,
        double &e,
        GBHypoParameters GBparameters
        )
{
    // normalized stress tensor
    double sigma_hat[NCOMP_TENSOR]; 
    // normalized deviatoric stress tenstor
    double sigma_hat_star[NCOMP_TENSOR];
    // Jaumann stress rate
    double sigma_o[NCOMP_TENSOR];
    // stress rate tensor
    double sigma_dot[NCOMP_TENSOR];
    // identity matrix upper
    double Identity[NCOMP_TENSOR] = {1.0,0,0,1.0,0,1.0};
    // void ratio
    //double e;

    //get tensors for computation of stress rate
    double tensor_Jaumann[AMREX_SPACEDIM*AMREX_SPACEDIM];
    double tensor_spin[AMREX_SPACEDIM*AMREX_SPACEDIM];
    double tensor_stress[AMREX_SPACEDIM*AMREX_SPACEDIM];

    // intermediate products
    double tensor_omega_dot_sigma[AMREX_SPACEDIM*AMREX_SPACEDIM];
    double tensor_sigma_dot_omega[AMREX_SPACEDIM*AMREX_SPACEDIM];

    // Linear contribution
    double L[NCOMP_TENSOR];
    // Nonlinear contribution
    double N[NCOMP_TENSOR];

    // GH hypoplastic model parameters
    double phi_c = GBparameters.phi_c; // critical state internal friction angle
    double h_s = GBparameters.h_s; // granular hardness
    double n = GBparameters.n; // pressure sensitive exponent : G paper says range (0.3, 0.5)
    double e_d0 = GBparameters.e_d0; // material constant for minimum void ratio
    double e_c0 = GBparameters.e_c0; // material constant for critical void ratio
    double e_i0 = GBparameters.e_i0; // material constant for maximum void ratio
    double alpha = GBparameters.alpha; //material constant (0.1, 0.3)  
    double beta = GBparameters.beta; //material constant (1.0, 1.1)
    double stress_cutoff = GBparameters.stress_cutoff;

    double trace_sigma = sigma[XX] + sigma[YY] + sigma[ZZ];
    
    if(trace_sigma>0){
        trace_sigma = 0;
    }

    if(trace_sigma != 0.0){
        // compute some varaibles
        double trace_eps_dot=eps_dot[XX]+eps_dot[YY]+eps_dot[ZZ];
        double modulus_eps_dot = sqrt(pow(trace_eps_dot,2));

        for(int i = 0; i < NCOMP_TENSOR; i++) {
            sigma_hat[i] = sigma[i]/trace_sigma;
            sigma_hat_star[i] = sigma_hat[i] - (1.0/3.0)*Identity[i];
        }
        double trace_sigma_hat_star = sigma_hat_star[XX]+sigma_hat_star[YY]+sigma_hat_star[ZZ];
        //double modulus_sigma = sqrt(pow(trace_sigma,2));
            

        double modulus_sigma_hat_star = sqrt(pow(trace_sigma_hat_star,2));

        double trace_sigma_hat_star_square = pow(sigma_hat_star[XX],2)+pow(sigma_hat_star[YY],2)+pow(sigma_hat_star[ZZ],2)
                                                    +2*pow(sigma_hat_star[XY],2)+2*pow(sigma_hat_star[YZ],2)+2*pow(sigma_hat_star[XZ],2);
        double det_sigma_hat_star = sigma_hat_star[XX]*sigma_hat_star[YY]*sigma_hat_star[ZZ]
                                        +2*sigma_hat_star[XY]*sigma_hat_star[YZ]*sigma_hat_star[XZ]
                                        -sigma_hat_star[XY]*sigma_hat_star[XY]*sigma_hat_star[ZZ]
                                        -sigma_hat_star[YZ]*sigma_hat_star[YZ]*sigma_hat_star[XX]
                                        -sigma_hat_star[XZ]*sigma_hat_star[XZ]*sigma_hat_star[YY];
        double trace_sigma_hat_star_cube = 3*det_sigma_hat_star;

        double trace_sigma_hat_dot_eps_dot = sigma_hat[XX]*eps_dot[XX] + sigma_hat[YY]*eps_dot[YY] +sigma_hat[ZZ]*eps_dot[ZZ]
                                            + 2*sigma_hat[XY]*eps_dot[XY] + 2*sigma_hat[YZ]*eps_dot[YZ] + 2*sigma_hat[XZ]*eps_dot[XZ];
        
        // compute c1 c2
        double c_1 = sqrt(3.0/8.0)*((3-sin(phi_c))/sin(phi_c));
        double c_2 = (3.0/8.0)*((3+sin(phi_c))/sin(phi_c));

        // update stiffness factor
        double h_i = 1/pow(c_1,2) + (1.0/3.0)*pow((e_i0+e_d0)/(e_c0-e_d0),alpha)*(1/(sqrt(3)*c_1));

        // update void ratio related variables
        double e_i = e_i0*exp(-pow(-trace_sigma/h_s,n)); // max void ratio
        double e_d = e_d0*exp(-pow(-trace_sigma/h_s,n)); // min void ratio
        double e_c = e_c0*exp(-pow(-trace_sigma/h_s,n)); // critical void ratio
        if(e>e_i){
            e = e_i;
        }
        else if (e<e_d){
            e = e_d;
        }

        // Pycnotropy factors
        double f_d = pow((e-e_d)/(e_c-e_d),alpha);
        double f_e = pow(e_c/e, beta);
        double f_b = pow(e_i0/e_c0, beta)*(h_s/(n*h_i))*((1+e_i)/e_i)*pow(-trace_sigma/h_s,1-n);
        double f_s = f_b*f_e;

        // compute lode angle
        double cos3theta = -sqrt(6.0)*trace_sigma_hat_star_cube/pow(trace_sigma_hat_star_square,1.5);

        double a_1 = 1.0/(c_1+c_2*modulus_sigma_hat_star*(1+cos3theta));

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
        //std::cout<< "sigma_dot: ["<< sigma_dot[XX]<<","<<sigma_dot[XY]<<","<<sigma_dot[XZ]<<","<<sigma_dot[YY]<<","<<sigma_dot[YZ]<<","<<sigma_dot[ZZ]<<"\n";
        //limit stress to be compressive
        sigma[XX] = (sigma[XX]<=stress_cutoff)?sigma[XX]:stress_cutoff;
        //sigma[YY] = (sigma[YY]<=stress_cutoff)?sigma[YY]:-0.0;
        //sigma[ZZ] = (sigma[ZZ]<=stress_cutoff)?sigma[ZZ]:-0.0;
    }
    else{
        if(e>e_i0){
            e = e_i0;
        }
        else if (e<e_d0){
            e = e_d0;
        }
    }
}