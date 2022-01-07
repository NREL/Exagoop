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
    void deposit_onto_grid(MultiFab& nodaldata,amrex::Array<amrex::Real,AMREX_SPACEDIM> gravity,int update_massvel,int update_forces);
    void interpolate_from_grid(MultiFab& nodaldata,
            int update_vel,int update_strainrate);

    void apply_constitutive_model(const amrex::Real& dt,
        amrex::Real E,amrex::Real v);

private:

};
#endif
