#ifndef CONSTITUTIVE_MODELS_H_
#define CONSTITUTIVE_MODELS_H_

#include<constants.H>
#include <AMReX.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void linear_elastic(amrex::Real eps[NCOMP_TENSOR],
        amrex::Real epsdot[NCOMP_TENSOR],
        amrex::Real sigma[NCOMP_TENSOR],
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Newtonian_Fluid(amrex::Real epsdot[NCOMP_TENSOR],
        amrex::Real sigma[NCOMP_TENSOR],
        amrex::Real dyn_visc,
		amrex::Real pressure)
{

    amrex::Real trace_epsdot=epsdot[XX]+epsdot[YY]+epsdot[ZZ];
    amrex::Real one_by_three=1.0/3.0;

    sigma[XX]=2.0*dyn_visc*(epsdot[XX]-one_by_three*trace_epsdot)-pressure;
    sigma[YY]=2.0*dyn_visc*(epsdot[YY]-one_by_three*trace_epsdot)-pressure;
    sigma[ZZ]=2.0*dyn_visc*(epsdot[ZZ]-one_by_three*trace_epsdot)-pressure;

    sigma[XY]=2.0*dyn_visc*(epsdot[XY]);
    sigma[XZ]=2.0*dyn_visc*(epsdot[XZ]);
    sigma[YZ]=2.0*dyn_visc*(epsdot[YZ]);
    for(int i=0; i<NCOMP_TENSOR; i++){
        amrex::Print()<<"\nsigma["<<i<<"] = "<<sigma[i];
    }

}

/**
 * @brief Yudong: implementation of Gudehus-Bauer hypoplastic model 
 * 
 * @param epsdot 
 * @param sigma 
 * @param e : void ratio for material properties
 * @param dt : time step for updating stress tensor
 * @return AMREX_GPU_DEVICE 
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void GB_hypoplastic(amrex::Real eps_dot[NCOMP_TENSOR],
        amrex::Real sigma[NCOMP_TENSOR],
        amrex::Real dt,
        amrex::Real &e)
{
    // normalized stress tensor
    amrex::Real sigma_hat[NCOMP_TENSOR]; 
    // normalized deviatoric stress tenstor
    amrex::Real sigma_hat_star[NCOMP_TENSOR];
    // stress rate tensor
    amrex::Real sigma_dot[NCOMP_TENSOR];
    // identity matrix upper
    amrex::Real Identity[NCOMP_TENSOR] = {1.0,0,0,1.0,0,1.0};
    // void ratio
    //amrex::Real e;

    // Linear contribution
    amrex::Real L[NCOMP_TENSOR];
    // Nonlinear contribution
    amrex::Real N[NCOMP_TENSOR];

    // hard code GB parameters for now. Should be inited from input files.
    const double phi_c = 30.0; // critical state internal friction angle
    const double h_s = 5800.0; // granular hardness
    const double n = 0.28; // pressure sensitive exponent : G paper says range (0.3, 0.5)
    const double e_d0 = 0.53; // material constant for minimum void ratio
    const double e_c0 = 0.84; // material constant for critical void ratio
    const double e_i0 = 1.0; // material constant for maximum void ratio
    const double alpha = 0.13; //material constant (0.1, 0.3)  
    const double beta = 1.0; //material constant (1.0, 1.1)

    // compute some varaibles
    amrex::Real trace_eps_dot=eps_dot[XX]+eps_dot[YY]+eps_dot[ZZ];
    //if(trace_eps_dot == 0.0){
    //    trace_eps_dot = TINYVAL; //avoid divide by 0
    //}
    amrex::Real modulus_eps_dot = sqrt(pow(trace_eps_dot,2));

    amrex::Real trace_sigma = sigma[XX] + sigma[YY] + sigma[ZZ];
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
    amrex::Real c_1 = sqrt(3.0/8.0)*((3-sin(phi_c)/sin(phi_c)));
    amrex::Real c_2 = (3.0/8.0)*((3+sin(phi_c)/sin(phi_c)));

    // update stiffness factor
    amrex::Real h_i = 1/pow(c_1,2) + (1.0/3.0)*pow((e_i0+e_d0)/(e_c0-e_d0),alpha)*(1/(sqrt(3)*c_1));

    if(trace_sigma>0){
        trace_sigma = 0;
    }

    if(trace_sigma != 0.0){
        // update void ratio related variables
        amrex::Real e_i = e_i0*exp(-pow(-trace_sigma/h_s,n)); // max void ratio
        amrex::Real e_d = e_d0*exp(-pow(-trace_sigma/h_s,n)); // min void ratio
        amrex::Real e_c = e_c0*exp(-pow(-trace_sigma/h_s,n)); // critical void ratio
        if(e>e_i){
            //amrex::Print()<<"\n Void ratio exeeds upper bond: "<<e<<"(e_i="<<e_i<<");";
            e = e_i;
        }
        else if (e<e_d){
            //amrex::Print()<<"\n Void ratio exeeds lower bond: "<<e<<"(e_d="<<e_d<<");";
            e = e_d;
        }

        // Pycnotropy factors
        amrex::Real f_d = pow((e-e_d)/(e_c-e_d),alpha);
        amrex::Real f_e = pow(e_c/e, beta);
        amrex::Real f_b = pow(e_i0/e_c0, beta)*(h_s/(n*h_i))*((1+e_i)/e_i)*pow(-trace_sigma/h_s,1-n);
        amrex::Real f_s = f_b*f_e;

        // compute lode angle
        amrex::Real cos3theta = -sqrt(6.0)*trace_sigma_hat_star_cube/pow(trace_sigma_hat_star_square,1.5);

        amrex::Real a_1 = 1.0/(c_1+c_2*modulus_sigma_hat_star*(1+cos3theta));
        /*
        for(int i=0; i<NCOMP_TENSOR; i++){
            amrex::Print()<<"\n____sigma_old["<<i<<"] = "<<sigma[i];
        }
        amrex::Print()<<"\ne = "<<e;
        amrex::Print()<<"\ntrace_sigma = "<<trace_sigma;
        amrex::Print()<<"\npower exponent for e_d = "<<-trace_sigma/h_s;
        amrex::Print()<<"\nexponent for e_d = "<<-pow(-trace_sigma/h_s,n);
        amrex::Print()<<"\ne_d = "<<e_d;
        amrex::Print()<<"\nf_b = "<<f_b;
        amrex::Print()<<"\nf_s = "<<f_s;
        amrex::Print()<<"\nf_d = "<<f_d;
        amrex::Print()<<"\nmodulus_eps_dot = "<<modulus_eps_dot;
        */
        // update stress based on stress rate
        for(int i=0; i<NCOMP_TENSOR; i++){
            N[i] = a_1*(sigma_hat[i]+sigma_hat_star[i]);
            // linear contribution
            L[i] = pow(a_1,2)*eps_dot[i]+sigma_hat[i]*trace_sigma_hat_dot_eps_dot;
            sigma_dot[i] = f_s*(L[i]+f_d*N[i]*modulus_eps_dot);
            sigma[i] = sigma[i]+sigma_dot[i]*dt;
            /*
            amrex::Print()<<"\nN["<<i<<"] = "<<N[i];
            amrex::Print()<<"\nL["<<i<<"] = "<<L[i];
            amrex::Print()<<"\nsigma_dot["<<i<<"] = "<<sigma_dot[i];
            amrex::Print()<<"\nsigma["<<i<<"] = "<<sigma[i];
            */
        }
    }
    else {
        //amrex::Print()<<"\n trace equals "<<trace_sigma;
        amrex::Print()<<"\n trace greater than or equals to zero , then set to : "<<trace_sigma;
        amrex::Real a_1 = 1.0/c_1;
        if(e>e_i0){
            //amrex::Print()<<"\n Void ratio exeeds upper bond: "<<e<<"(e_i0="<<e_i0<<");";
            e = e_i0;
        }
        else if (e<e_d0){
            //amrex::Print()<<"\n Void ratio exeeds lower bond: "<<e<<"(e_d0="<<e_d0<<");";
            e = e_d0;
        }
        amrex::Real f_e = pow(e_c0/e, beta);// modified e_c to equal to e_c0 when trace_sigma == 0;
        amrex::Real f_b = pow(e_i0/e_c0, beta)*(h_s/(n*h_i));
        amrex::Real f_s = f_b*f_e;

        for(int i=0; i<NCOMP_TENSOR; i++){
            L[i] = pow(a_1,2)*eps_dot[i];
            sigma_dot[i] = f_s*L[i];
            sigma[i] = sigma[i]+sigma_dot[i]*dt;
            //amrex::Print()<<"\nsigma["<<i<<"] = "<<sigma[i];
        }
    }

}
#endif
