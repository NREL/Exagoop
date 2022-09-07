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

    void InitParticles (const std::string & filename,Real &total_mass, Real &total_vol);

    void InitParticles (Real mincoords[AMREX_SPACEDIM],Real maxcoords[AMREX_SPACEDIM], 
        Real vel[AMREX_SPACEDIM],
        Real dens, int constmodel, 
        Real E, Real nu,Real bulkMod, Real Gama_pres,Real visc,
        int do_multi_part_per_cell,Real &total_mass,Real &total_vol, Real inital_void_ratio);

    void moveParticles (const amrex::Real& dt,
                        int bclo[AMREX_SPACEDIM],
                        int bchi[AMREX_SPACEDIM],
                        int lsetbc,
                        amrex::Real wall_mu_lo[AMREX_SPACEDIM],
                        amrex::Real wall_mu_hi[AMREX_SPACEDIM],
                        amrex::Real wall_vel_lo[AMREX_SPACEDIM*AMREX_SPACEDIM],
                        amrex::Real wall_vel_hi[AMREX_SPACEDIM*AMREX_SPACEDIM],
                        amrex::Real lset_wall_mu);

    void updateVolume (const amrex::Real& dt);
    amrex::Real CalculateExactVelocity(int modenumber,amrex::Real E, amrex::Real rho, amrex::Real v0,amrex::Real L, amrex::Real time);
    void writeParticles (std::string prefix_particlefilename, int num_of_digits_in_filenames, const int n);
    void writeCheckpointFile(std::string prefix_particlefilename, int num_of_digits_in_filenames, amrex::Real cur_time,  int nstep, int output_it);
    void WriteHeader(const std::string& name, bool is_checkpoint, amrex::Real cur_time, int nstep, int EB_generate_max_level, int output_it) const;
    void readCheckpointFile(std::string & restart_chkfile, int &nstep, double &cur_time, int &output_it);
    void deposit_onto_grid(MultiFab& nodaldata,
                           Array<Real,AMREX_SPACEDIM> gravity,
                           int external_loads_present,
                           Array<Real,AMREX_SPACEDIM> force_slab_lo,
                           Array<Real,AMREX_SPACEDIM> force_slab_hi,
                           Array<Real,AMREX_SPACEDIM> extforce,
                           int update_massvel,int update_forces, 
                           amrex::Real mass_tolerance,
						   Array<int, AMREX_SPACEDIM> order_scheme_directional,
						   Array<int, AMREX_SPACEDIM> periodic);

    void interpolate_from_grid(MultiFab& nodaldata,
                               int update_vel,int update_strainrate,
							   Array <int,AMREX_SPACEDIM> order_scheme_directional,
							   Array <int,AMREX_SPACEDIM> periodic,
							   amrex::Real alpha_pic_flip);

    void move_particles_from_nodevel(MultiFab& nodaldata,const amrex::Real& dt,
            int bclo[AMREX_SPACEDIM],int bchi[AMREX_SPACEDIM],
            int order_scheme);
    
    void CalculateEnergies(Real &TKE,Real &TSE);
    void CalculateVelocity(Real &Vcm);
    void CalculateVelocityCantilever(Real &Vcm);
    void WriteDeflectionCantilever();
    void FindWaterFront(Real &Xwf);

    void apply_constitutive_model(MPMspecs specs, const amrex::Real& dt,amrex::Real applied_strainrate);

    void interpolate_mass_from_grid(MultiFab& nodaldata,int order_scheme);

    void update_density_field(MultiFab& nodaldata,int refratio,amrex::Real smoothfactor);
    amrex::Real Calculate_time_step(amrex::Real CFL, amrex::Real dtmax,amrex::Real dtmin);
    void removeParticlesInsideEB();

private:

    ParticleType generate_particle
        (Real x,Real y,Real z,
        Real vel[AMREX_SPACEDIM],
        Real dens, Real vol, int constmodel, Real E, Real nu,
        Real bulkmod, Real Gama_pres,Real visc, Real initial_void_ratio);

};
#endif
