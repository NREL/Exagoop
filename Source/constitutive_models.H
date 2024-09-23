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
void linear_elastic(amrex::Real delta_eps[NCOMP_TENSOR],
                    amrex::Real delta_sigma[NCOMP_TENSOR],
                    amrex::Real E, amrex::Real v)
{
    amrex::Real const1=E/(one+v)/(one-two*v);
    amrex::Real const2=E/(one+v);

    delta_sigma[XX]=const1*((1-v)*delta_eps[XX] + v*(delta_eps[YY]+delta_eps[ZZ]));
    delta_sigma[YY]=const1*((1-v)*delta_eps[YY] + v*(delta_eps[XX]+delta_eps[ZZ]));
    delta_sigma[ZZ]=const1*((1-v)*delta_eps[ZZ] + v*(delta_eps[XX]+delta_eps[YY]));

    delta_sigma[XY]=const2*delta_eps[XY];
    delta_sigma[XZ]=const2*delta_eps[XZ];
    delta_sigma[YZ]=const2*delta_eps[YZ];
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

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Hyper_Elastic_Solid(amrex::Real deform_grad[NCOMP_TENSOR],
		amrex::Real sigma[NCOMP_TENSOR],
		amrex::Real jacob,
        amrex::Real mu,
        amrex::Real lambda
        )
{
	amrex::Real one_by_J = 1.0/jacob;
	amrex::Real lnJ = log(jacob);
	amrex::Real B[NCOMP_TENSOR];
	//F*Ftranspose
	B[XX] = deform_grad[XX]*deform_grad[XX]+deform_grad[XY]*deform_grad[XY]+deform_grad[XZ]*deform_grad[XZ];
	B[YY] = deform_grad[XY]*deform_grad[XY]+deform_grad[YY]*deform_grad[YY]+deform_grad[YZ]*deform_grad[YZ];
	B[ZZ] = deform_grad[XZ]*deform_grad[XZ]+deform_grad[YZ]*deform_grad[YZ]+deform_grad[ZZ]*deform_grad[ZZ];

	B[XY] = deform_grad[XX]*deform_grad[XY]+deform_grad[XY]*deform_grad[YY]+deform_grad[XZ]*deform_grad[YZ];
	B[XZ] = deform_grad[XX]*deform_grad[XZ]+deform_grad[XY]*deform_grad[YZ]+deform_grad[XZ]*deform_grad[ZZ];
	B[YZ] = deform_grad[XY]*deform_grad[XZ]+deform_grad[YY]*deform_grad[YZ]+deform_grad[YZ]*deform_grad[ZZ];

	sigma[XX] = one_by_J*(mu*(B[XX]-1.0)+lambda*lnJ);
	sigma[YY] = one_by_J*(mu*(B[YY]-1.0)+lambda*lnJ);
	sigma[ZZ] = one_by_J*(mu*(B[ZZ]-1.0)+lambda*lnJ);

	sigma[XY] = one_by_J*(mu*B[XY]);
	sigma[YZ] = one_by_J*(mu*B[YZ]);
	sigma[XZ] = one_by_J*(mu*B[XZ]);

}


#endif
