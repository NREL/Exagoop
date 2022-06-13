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

    void InitParticles (const std::string & filename,Real *total_mass, Real *total_vol);
    void moveParticles (const amrex::Real& dt);
    void updatevolume (const amrex::Real& dt, const amrex::Real& K, const amrex::Real& G );
    void writeParticles (const int n);  
    void deposit_onto_grid(MultiFab& nodaldata,
                           Array<Real,AMREX_SPACEDIM> gravity,
                           int external_loads_present,
                           Array<Real,AMREX_SPACEDIM> force_slab_lo,
                           Array<Real,AMREX_SPACEDIM> force_slab_hi,
                           Array<Real,AMREX_SPACEDIM> extforce,
                           int update_massvel,int update_forces, amrex::Real mass_tolerance);
    void interpolate_from_grid(MultiFab& nodaldata,
                               int update_vel,int update_strainrate,int order,amrex::Real alpha_pic_flip);
    void move_particles_from_nodevel(MultiFab& nodaldata,const amrex::Real& dt,int order_scheme);

    void apply_constitutive_model(const amrex::Real& dt,
                                  amrex::Real E,
								  amrex::Real v,
								  int Constitutive_Model,
								  amrex::Real K,
								  amrex::Real G,
								  amrex::Real dyn_visc,
								  amrex::Real applied_strainrate
								  );
    void interpolate_mass_from_grid(MultiFab& nodaldata,int order_scheme);

    void update_density_field(MultiFab& nodaldata,int refratio,amrex::Real smoothfactor);
    amrex::Real Calculate_time_step(amrex::Real K);

private:

};
#endif
