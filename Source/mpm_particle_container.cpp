#include <mpm_particle_container.H>
#include <interpolants.H>
#include <constitutive_models.H>
#include <gb_hypoplastic.H>
#include <miuIrheology.H>

using namespace amrex;

void MPMParticleContainer::apply_constitutive_model_on_rigidparticles(const amrex::Real& dt,
                                                    amrex::Real applied_strainrate=0.0)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);
        Real pinf=0.0;

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        int np = aos.numRealParticles();
        int ng = aos.numNeighborParticles();
        int nt = np+ng;

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(nt,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            Real p_inf=0.0;
            ParticleType& p = pstruct[i];

            if(p.idata(intData::phase)==1)
            {
                amrex::Real xp[AMREX_SPACEDIM];
                amrex::Real strainrate[NCOMP_TENSOR];
                amrex::Real strain[NCOMP_TENSOR];
                amrex::Real stress[NCOMP_TENSOR];

                for(int d=0;d<NCOMP_TENSOR;d++)
                {
                    p.rdata(realData::strain+d) += dt*p.rdata(realData::strainrate+d);
                }
                //apply axial strain
                p.rdata(realData::strain+XX) += dt*applied_strainrate;
                p.rdata(realData::strain+YY) += dt*applied_strainrate;
                p.rdata(realData::strain+ZZ) += dt*applied_strainrate;

                for(int d=0;d<NCOMP_TENSOR;d++)
                {
                    strainrate[d]=p.rdata(realData::strainrate+d);
                    strain[d]=p.rdata(realData::strain+d);
                }

                if(p.idata(intData::constitutive_model)==0)             //Elastic solid
                {
                    linear_elastic(strain,strainrate,stress,p.rdata(realData::E),p.rdata(realData::nu));
                }
                else if(p.idata(intData::constitutive_model)==1)                //Viscous fluid with approximate EoS
                {
                    p.rdata(realData::pressure) = p.rdata(realData::Bulk_modulus)*
                    (pow(1/p.rdata(realData::jacobian),p.rdata(realData::Gama_pressure))-1.0)+p_inf;
                    Newtonian_Fluid(strainrate,stress,p.rdata(realData::Dynamic_viscosity),p.rdata(realData::pressure));
                }

                for(int d=0;d<NCOMP_TENSOR;d++)
                {
                    p.rdata(realData::stress+d)=stress[d];
                }
            }
        });
    }
}

void MPMParticleContainer::apply_constitutive_model(const amrex::Real& dt,
                                                    amrex::Real applied_strainrate=0.0)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);
        Real pinf=0.0;

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        int np = aos.numRealParticles();
        int ng = aos.numNeighborParticles();
        int nt = np+ng;

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(nt,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            Real p_inf=0.0;
            ParticleType& p = pstruct[i];

            if(p.idata(intData::phase)==0)
            {
                amrex::Real xp[AMREX_SPACEDIM];
                amrex::Real strainrate[NCOMP_TENSOR];
                amrex::Real strain[NCOMP_TENSOR];
                amrex::Real stress[NCOMP_TENSOR];

                for(int d=0;d<NCOMP_TENSOR;d++)
                {
                    p.rdata(realData::strain+d) += dt*p.rdata(realData::strainrate+d);
                }
                //apply axial strain
                p.rdata(realData::strain+XX) += dt*applied_strainrate;
                p.rdata(realData::strain+YY) += dt*applied_strainrate;
                p.rdata(realData::strain+ZZ) += dt*applied_strainrate;

                for(int d=0;d<NCOMP_TENSOR;d++)
                {
                    strainrate[d]=p.rdata(realData::strainrate+d);
                    strain[d]=p.rdata(realData::strain+d);
                }

                if(p.idata(intData::constitutive_model)==0)		//Elastic solid
                {
                    linear_elastic(strain,strainrate,stress,p.rdata(realData::E),p.rdata(realData::nu));
                }
                else if(p.idata(intData::constitutive_model)==1)		//Viscous fluid with approximate EoS
                {
                    p.rdata(realData::pressure) = p.rdata(realData::Bulk_modulus)*
                    (pow(1/p.rdata(realData::jacobian),p.rdata(realData::Gama_pressure))-1.0)+p_inf;
                    Newtonian_Fluid(strainrate,stress,p.rdata(realData::Dynamic_viscosity),p.rdata(realData::pressure));
                }

                for(int d=0;d<NCOMP_TENSOR;d++)
                {
                    p.rdata(realData::stress+d)=stress[d];
                }
            }
        });
    }
} 

void MPMParticleContainer::apply_constitutive_model_delta(const amrex::Real& dt,
                                                          amrex::Real applied_strainrate=0.0)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);
        Real pinf=0.0;

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        int np = aos.numRealParticles();
        int ng = aos.numNeighborParticles();
        int nt = np+ng;

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(nt,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            Real p_inf=0.0;
            ParticleType& p = pstruct[i];

            if(p.idata(intData::phase)==0)
            {
				amrex::Real xp[AMREX_SPACEDIM];
				amrex::Real strainrate[NCOMP_TENSOR];
				amrex::Real delta_strain[NCOMP_TENSOR];
				amrex::Real delta_stress[NCOMP_TENSOR];

				for(int d=0;d<NCOMP_TENSOR;d++)
				{
					p.rdata(realData::strain+d) += dt*p.rdata(realData::strainrate+d);
				}
				//apply axial strain
				p.rdata(realData::strain+XX) += dt*applied_strainrate;
				p.rdata(realData::strain+YY) += dt*applied_strainrate;
				p.rdata(realData::strain+ZZ) += dt*applied_strainrate;

				for(int d=0;d<NCOMP_TENSOR;d++)
				{
					delta_strain[d]=dt*p.rdata(realData::strainrate+d);
				}

				//apply axial strain
				delta_strain[XX] += dt*applied_strainrate;
				delta_strain[YY] += dt*applied_strainrate;
				delta_strain[ZZ] += dt*applied_strainrate;

				if(p.idata(intData::constitutive_model)==0)		//Elastic solid
				{
					linear_elastic(delta_strain,delta_stress,p.rdata(realData::E),p.rdata(realData::nu));
				}
				else if(p.idata(intData::constitutive_model)==1)		//Viscous fluid with approximate EoS
				{
					amrex::Abort("\nDelta strain model for weakly compressible fluids not implemented yet.");
				}

				for(int d=0;d<NCOMP_TENSOR;d++)
				{
					p.rdata(realData::stress+d)+=delta_stress[d];
				}
            }
        });
    }
}

//yli add for gbhypo
void MPMParticleContainer::apply_constitutive_model(MPMspecs specs,
													const amrex::Real& dt,
                                                    amrex::Real applied_strainrate=0.0)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);
        Real pinf=0.0;

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        int np = aos.numRealParticles();
        int ng = aos.numNeighborParticles();
        int nt = np+ng;

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(nt,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            Real p_inf=0.0;
            ParticleType& p = pstruct[i];

            if(p.idata(intData::phase)==0)
            {
				amrex::Real xp[AMREX_SPACEDIM];
				amrex::Real strainrate[NCOMP_TENSOR];
				//yli add for gbhypo
				amrex::Real spinrate[NCOMP_TENSOR];
				//end yli add
				amrex::Real strain[NCOMP_TENSOR];
				amrex::Real stress[NCOMP_TENSOR];

				for(int d=0;d<NCOMP_TENSOR;d++)
				{
					p.rdata(realData::strain+d) += dt*p.rdata(realData::strainrate+d);
				}
				//apply axial strain
				p.rdata(realData::strain+XX) += dt*applied_strainrate;
				p.rdata(realData::strain+YY) += dt*applied_strainrate;
				p.rdata(realData::strain+ZZ) += dt*applied_strainrate;

				for(int d=0;d<NCOMP_TENSOR;d++)
				{
					strainrate[d]=p.rdata(realData::strainrate+d);
					strain[d]=p.rdata(realData::strain+d);
					//yli add for gbhypo
					spinrate[d]=p.rdata(realData::spinrate+d);
					stress[d]=p.rdata(realData::stress+d);
					//end yli add

				}

				if(p.idata(intData::constitutive_model)==0)		//Elastic solid
				{
					linear_elastic(strain,strainrate,stress,p.rdata(realData::E),p.rdata(realData::nu));
				}
				else if(p.idata(intData::constitutive_model)==1)		//Viscous fluid with approximate EoS
				{
					p.rdata(realData::pressure) = p.rdata(realData::Bulk_modulus)*
						(pow(1/p.rdata(realData::jacobian),p.rdata(realData::Gama_pressure))-1.0)+p_inf;
					Newtonian_Fluid(strainrate,stress,p.rdata(realData::Dynamic_viscosity),p.rdata(realData::pressure));
				}
				//yli add for gbhypo
				else if(p.idata(intData::constitutive_model)==2)		//Yudong: GB hypoplastic model for granular flow
            	{
					//GB hypoplastic model here.
					//amrex::Real e = p.rdata(realData::void_ratio);
					//amrex::Real e_before = e;

					//GB_hypoplastic(strainrate,spinrate,stress,dt, e,specs);
                    GB_hypoplastic(strainrate,spinrate,stress,dt, p,specs);
					//p.rdata(realData::void_ratio) = e;
            	}
                else if(p.idata(intData::constitutive_model)==3) //miu-I rheology
                {
                    miuI(strainrate,stress,p,specs);
                }
				//end yli add

				for(int d=0;d<NCOMP_TENSOR;d++)
				{
					p.rdata(realData::stress+d)=stress[d];
				}
            }
        });
    }
} 

//end yli add