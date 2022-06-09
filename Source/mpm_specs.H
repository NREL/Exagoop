#ifndef MPM_PARTICLE_SPECS_H_
#define MPM_PARTICLE_SPECS_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

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
	//vol_init=27,
    count=27
  };
};

struct intData
{
  enum 
  {
    phase = 0,
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

        int Constitutive_Model = 1;		//0-->linear elastic, 1-->NEwtonian fluid
        int order_scheme=3;				//1-->tent function, 3--> cubic spline shape function
        Real Dynamic_Viscosity=0.001;	//In case fluid simulations are involved
        Real Bulk_Modulous = 2e6;		//In case fluid simulations are involved
        Real Youngs_modulus=1e9;
        Real Poissons_ratio=0.1;
        Real applied_strainrate=0.0;
        Real applied_strainrate_time=0.0;
        int dens_field_output=0;
        int dens_field_gridratio=1;
        Real smoothfactor=1.0;
        int external_loads_present=0.0;
        
        void read_mpm_specs()
        {
            ParmParse pp("mpm");
            pp.get("prob_lo",plo);
            pp.get("prob_hi",phi);
            pp.get("ncells", ncells);
            pp.get("is_it_periodic",periodic);
            pp.get("max_grid_size", max_grid_size);
         
            pp.get("max_steps", maxsteps);
            pp.get("final_time",final_time);
            pp.get("timestep", timestep);

            pp.query("gravity", gravity);

            write_output_time=final_time/10.0;
            pp.query("write_output_time", write_output_time);

            screen_output_time=write_output_time/10.0;
            pp.query("screen_output_time", screen_output_time);

            pp.query("num_redist", num_redist);

            particlefilename="mpm_particles.dat";
            pp.query("particle_file",particlefilename);

            pp.query("Youngs_modulus",Youngs_modulus);
            pp.query("Poissons_ratio",Poissons_ratio);
            pp.query("Constitutive_Model",Constitutive_Model);
            amrex::Print()<<"\nCons model = "<<Constitutive_Model;

            pp.query("Dynamic_Viscosity",Dynamic_Viscosity);
            pp.query("Bulk_Modulous",Bulk_Modulous);

            pp.query("applied_strainrate",applied_strainrate);
            applied_strainrate_time=final_time;
            pp.query("applied_strainrate_time",applied_strainrate_time);

            pp.query("order_scheme",order_scheme);

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
        }
};
#endif
