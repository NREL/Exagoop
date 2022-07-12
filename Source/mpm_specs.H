#ifndef MPM_PARTICLE_SPECS_H_
#define MPM_PARTICLE_SPECS_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <constants.H>

using namespace amrex;

struct realData
{
    enum 
    { // Particle data stored in p.m_data
        radius=0,
        xvel=1,
        yvel=2,
        zvel=3,
        strainrate=4,
        strain=10,
        stress=16,
        volume=22,
        mass=23,
        density=24,
        jacobian=25,
        pressure=26,
        vol_init=27,
        E=28,
        nu=29,
        Bulk_modulous=30,
        Gama_pressure=31,
        Dynamic_viscosity=32,
        void_ratio=33, // for hypoplastic model
        count
    };
};

struct intData
{
    enum 
    {
        phase = 0,
        constitutive_model,
        count
    };
};

class MPMspecs
{
    public:

        Array<Real,AMREX_SPACEDIM> plo{AMREX_D_DECL(0.0, 0.0, 0.0)};
        Array<Real,AMREX_SPACEDIM> phi{AMREX_D_DECL(0.0, 0.0, 0.0)};
        Array<int,AMREX_SPACEDIM> ncells{AMREX_D_DECL(0, 0, 0)};
        Array<int,AMREX_SPACEDIM> periodic{AMREX_D_DECL(0, 0, 0)};
        Array<Real,AMREX_SPACEDIM> gravity{AMREX_D_DECL(0.0, -9.81, 0.0)};

        Array<Real,AMREX_SPACEDIM> force_slab_lo{AMREX_D_DECL(0.0, 0.0, 0.0)};
        Array<Real,AMREX_SPACEDIM> force_slab_hi{AMREX_D_DECL(0.0, 0.0, 0.0)};
        Array<Real,AMREX_SPACEDIM> extforce{AMREX_D_DECL(0.0, 0.0, 0.0)};

        Real timestep=0.0;
        int num_redist=10;
        int max_grid_size=8;
        Real final_time;
        Real write_output_time;
        Real screen_output_time;
        int maxsteps=1;
        std::string particlefilename;

        //1-->tent function, 3--> cubic spline shape function
        int order_scheme=1;				

        Real applied_strainrate=0.0;
        Real applied_strainrate_time=0.0;
        int dens_field_output=0;
        int dens_field_gridratio=1;
        Real smoothfactor=1.0;
        int external_loads_present=0.0;
        int fixed_timestep = 1;
        Real dtmin=1e-5;
        Real CFL=0.1;
        Real mass_tolerance = 1e-8;
        Real alpha_pic_flip = 0.95;
        Real total_mass = 0.0;
        Real total_vol = 0.0;
        int stress_update_scheme=1;
        int print_diagnostics=0;
        
        Vector<int> bclo;
        Vector<int> bchi;

        void read_mpm_specs()
        {
            bclo.resize(AMREX_SPACEDIM);
            bchi.resize(AMREX_SPACEDIM);
            

            ParmParse pp("mpm");
            pp.get("prob_lo",plo);
            pp.get("prob_hi",phi);
            pp.get("ncells", ncells);
            pp.get("is_it_periodic",periodic);
            pp.get("max_grid_size", max_grid_size);

            pp.get("max_steps", maxsteps);
            pp.get("final_time",final_time);

            pp.query("gravity", gravity);

            write_output_time=final_time/10.0;
            pp.query("write_output_time", write_output_time);

            screen_output_time=write_output_time/10.0;
            pp.query("screen_output_time", screen_output_time);

            pp.query("num_redist", num_redist);

            particlefilename="mpm_particles.dat";
            pp.query("particle_file",particlefilename);

            pp.query("fixed_timestep",fixed_timestep);
            pp.query("stress_update_scheme",stress_update_scheme);
            if(fixed_timestep==1)
            {
                pp.get("timestep", timestep);
            }
            else
            {
                timestep =0.0;
            }

            pp.query("alpha_pic_flip",alpha_pic_flip);
            pp.query("dtmin",dtmin);
            pp.query("CFL",CFL);

            pp.query("applied_strainrate",applied_strainrate);
            applied_strainrate_time=final_time;
            pp.query("applied_strainrate_time",applied_strainrate_time);

            pp.query("order_scheme",order_scheme);
            pp.query("mass_tolerance",mass_tolerance);

            pp.query("dens_field_output",dens_field_output);
            pp.query("smoothfactor",smoothfactor);
            if(dens_field_output==1)
            {
                pp.get("dens_field_gridratio",dens_field_gridratio);
            }

            pp.query("external_loads",external_loads_present); 
            if(external_loads_present)
            {  
                pp.get("force_slab_lo", force_slab_lo);
                pp.get("force_slab_hi", force_slab_hi);
                pp.get("extforce", extforce);
            }
            pp.query("print_diagnostics",print_diagnostics);
            
            //by default it is periodic 
            bclo[0]=BC_SLIPWALL;
            bclo[1]=BC_SLIPWALL;
            bclo[2]=BC_SLIPWALL;

            bchi[0]=BC_SLIPWALL;
            bchi[1]=BC_SLIPWALL;
            bchi[2]=BC_SLIPWALL;
            
            pp.queryarr("bc_lower",bclo);
            pp.queryarr("bc_upper",bchi);
        }
};
#endif
