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
        Bulk_modulus=30,
        Gama_pressure=31,
        Dynamic_viscosity=32,
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
        std::string restart_checkfile;

        //Output parameters
        std::string prefix_particlefilename = "plt";
        std::string prefix_gridfilename =  "nplt";
        std::string prefix_densityfilename = "dens";
        std::string prefix_checkpointfilename = "chk";
        int num_of_digits_in_filenames=6;

        //Diagnostic parameters
        int print_diagnostics=0;
        bool is_standard_test = true;
        int test_number = 0;	//0->is_standard_test=false,1-->axial vibration of bar, 2-->Dam break, 3--> Elastic collision of disks, 4--> Lateral deflection of beam

        	//Parameters for Axial bar test
        	amrex::Real axial_bar_E = 100.0;
        	amrex::Real axial_bar_rho = 1.0;
        	amrex::Real axial_bar_L = 25.0;
        	int axial_bar_modenumber = 1;
        	amrex::Real axial_bar_v0 = 0.1;

        	//Parameters for dam break
        	amrex::Real dam_break_H1 = 0.2;
        	amrex::Real dam_break_g = 9.81;


        //1-->tent function, 3--> cubic spline shape function
        int order_scheme=1;				

        Real applied_strainrate=0.0;
        Real applied_strainrate_time=0.0;
        int dens_field_output=0;
        int dens_field_gridratio=1;
        Real smoothfactor=1.0;
        int external_loads_present=0.0;
        int fixed_timestep = 0;				//Use fixed time step provided by user if this flag is 1
        Real dt_max_limit=1e0;
        Real dt_min_limit=1e-12;
        Real CFL=0.1;
        Real mass_tolerance = 1e-8;
        Real alpha_pic_flip = 0.95;
        Real total_mass = 0.0;
        Real total_vol = 0.0;
        int stress_update_scheme=1;
        

        Vector<int> bclo;
        Vector<int> bchi;

        Array<int,AMREX_SPACEDIM> order_scheme_directional{AMREX_D_DECL(1,1,1)};
       
        int use_autogen=0; 
        Vector<Real> autogen_mincoords; 
        Vector<Real> autogen_maxcoords; 
        Vector<Real> autogen_vel; 
        int autogen_multi_part_per_cell=1;
        Real autogen_dens=1000.0;
        int autogen_constmodel=1;
        Real autogen_E=1e9;
        Real autogen_nu=0.1;
        Real autogen_bulkmod=0.0;
        Real autogen_Gama_pres=1.4;
        Real autogen_visc=0.001;

        int levelset_bc=1; //no slip wall
        Real levelset_wall_mu=0.1;


        Array<Real,AMREX_SPACEDIM> wall_mu_lo{AMREX_D_DECL(0.1, 0.1, 0.1)};
        Array<Real,AMREX_SPACEDIM> wall_mu_hi{AMREX_D_DECL(0.1, 0.1, 0.1)};

        Array<Real,AMREX_SPACEDIM*AMREX_SPACEDIM> wall_vel_lo{0.0};
        Array<Real,AMREX_SPACEDIM*AMREX_SPACEDIM> wall_vel_hi{0.0};
        
        Real levset_mu=0.1;

        void read_mpm_specs()
        {
            bclo.resize(AMREX_SPACEDIM);
            bchi.resize(AMREX_SPACEDIM);

            
            autogen_mincoords.resize(AMREX_SPACEDIM);
            autogen_maxcoords.resize(AMREX_SPACEDIM);
            autogen_vel.resize(AMREX_SPACEDIM);
            
            ParmParse pp1("amr");
            restart_checkfile="";
            pp1.query("restart_checkfile",restart_checkfile);



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

            pp.query("use_autogen",use_autogen);

            //Read output parameters
            pp.query("prefix_particlefilename",prefix_particlefilename);
            pp.query("prefix_gridfilename",prefix_gridfilename);
            pp.query("prefix_densityfilename",prefix_densityfilename);
            pp.query("prefix_checkpointfilename",prefix_checkpointfilename);
            pp.query("num_of_digits_in_filenames",num_of_digits_in_filenames);

            //Reading diagnostic parameters
            pp.query("is_standard_test",is_standard_test);
            pp.query("test_number",test_number);
            pp.query("axial_bar_E",axial_bar_E);
            pp.query("axial_bar_rho",axial_bar_rho);
            pp.query("axial_bar_L",axial_bar_L);
            pp.query("axial_bar_modenumber",axial_bar_modenumber);
            pp.query("axial_bar_v0",axial_bar_v0);


            if(!use_autogen)
            {
                particlefilename="mpm_particles.dat";
                pp.query("particle_file",particlefilename);
            }
            else
            {
                pp.getarr("mincoords_autogen",autogen_mincoords);
                pp.getarr("maxcoords_autogen",autogen_maxcoords);
                pp.getarr("vel_autogen",autogen_vel);
                pp.query("constmodel_autogen",autogen_constmodel);
                pp.query("dens_autogen",autogen_dens);
                pp.query("E_autogen",autogen_constmodel);
                pp.query("nu_autogen",autogen_constmodel);
                pp.query("bulkmod_autogen",autogen_bulkmod);
                pp.query("Gama_pres_autogen",autogen_Gama_pres);
                pp.query("visc_autogen",autogen_visc);
                pp.query("multi_part_per_cell_autogen",autogen_multi_part_per_cell);
            }

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
            pp.query("dt_max_limit",dt_max_limit);
            pp.query("dt_min_limit",dt_min_limit);
            pp.query("CFL",CFL);

            pp.query("applied_strainrate",applied_strainrate);
            applied_strainrate_time=final_time;
            pp.query("applied_strainrate_time",applied_strainrate_time);

            pp.query("order_scheme",order_scheme);
            pp.query("mass_tolerance",mass_tolerance);

            pp.query("dens_field_output",dens_field_output);
            pp.query("smoothfactor",smoothfactor);
            pp.query("dens_field_gridratio",dens_field_gridratio);

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

            pp.query("levelset_bc",levelset_bc);
            pp.query("levelset_wall_mu",levelset_wall_mu);


            pp.query("wall_mu_lo",wall_mu_lo);
            pp.query("wall_mu_hi",wall_mu_hi);

            for(int i=0;i<AMREX_SPACEDIM*AMREX_SPACEDIM;i++)
            {
                wall_vel_lo[i]=0.0;
                wall_vel_hi[i]=0.0;
            }

            pp.query("wall_vel_lo",wall_vel_lo);
            pp.query("wall_vel_hi",wall_vel_hi);
        }
};
#endif
