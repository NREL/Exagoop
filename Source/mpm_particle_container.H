#ifndef MPM_PARTICLE_CONTAINER_H_
#define MPM_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>
#include <mpm_specs.H>
#include <constants.H>

class MPMParticleContainer
    : public amrex::NeighborParticleContainer<realData::count, intData::count>
{

public:

    MPMParticleContainer (const amrex::Geometry            geom,
                         const amrex::DistributionMapping   dmap,
                         const amrex::BoxArray              ba,
                         int                                numcells)

    : NeighborParticleContainer<realData::count, intData::count>(geom, dmap, ba, numcells)
    {
    
    }

    void InitParticles (const std::string & filename);
    void moveParticles (const amrex::Real& dt);
    void writeParticles (const int n);  
    void deposit_onto_grid(MultiFab& nodaldata,
                           Array<Real,AMREX_SPACEDIM> gravity,
                           int external_loads_present,
                           Array<Real,AMREX_SPACEDIM> force_slab_lo,
                           Array<Real,AMREX_SPACEDIM> force_slab_hi,
                           Array<Real,AMREX_SPACEDIM> extforce,
                           int update_massvel,int update_forces); 
    void interpolate_from_grid(MultiFab& nodaldata,
                               int update_vel,int update_strainrate);

    void apply_constitutive_model(const amrex::Real& dt,
                                  amrex::Real E,
								  amrex::Real v,
								  int Constitutive_Model,
								  amrex::Real dyn_visc,
								  amrex::Real applied_strainrate
								  );

    void update_density_field(MultiFab& nodaldata,int refratio,amrex::Real smoothfactor);

private:

};
#endif
