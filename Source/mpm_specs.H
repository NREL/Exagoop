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
        xvel_prime=4,                           //For position update. Check the review paper page 34
        yvel_prime=5,
        zvel_prime=6,
        strainrate=7, 				//4
        strain=13,	  			//10
        stress=19,				//16
        deformation_gradient=25,	        //22
        volume=34,				//31
        mass=35,				//32
        density=36,				//33
        jacobian=37,				//34
        pressure=38,				//35
        vol_init=39,				//36
        E=40,					//37
        nu=41,					//38
        Bulk_modulus=42,			//39
        Gama_pressure=43,			//40
        Dynamic_viscosity=44,		        //41
        yacceleration=45,			//45	This is for storing the acceleration of each particle in the y-direction. This variable is not generic and used only for the spring-mass and membrane simulations only
        count
    };
};

struct intData
{
    enum 
    {
        phase = 0,
        rigid_body_id,
        constitutive_model,
        count
    };
};

struct Rigid_Bodies
{
	int Rigid_Body_Id;
	int num_of_mp;
	amrex::Real total_mass;
	amrex::Real total_volume;
	Array <amrex::Real,AMREX_SPACEDIM> gravity;
	Array <amrex::Real,AMREX_SPACEDIM> force_external;
	Array <amrex::Real,AMREX_SPACEDIM> force_internal;
	Array <amrex::Real,AMREX_SPACEDIM> velocity;
	Array <amrex::Real,AMREX_SPACEDIM> imposed_velocity;
	amrex::Real Damping_Coefficient;
	//Body forces
	int position_update_method;			//0-->Constant Velocity 1-->Particle dynamics solver
	bool enable_weight;
	bool enable_damping_force;
};



class MPMspecs
{
    public:

        Array<Real,AMREX_SPACEDIM> plo{AMREX_D_DECL(0.0, 0.0, 0.0)};
        Array<Real,AMREX_SPACEDIM> phi{AMREX_D_DECL(0.0, 0.0, 0.0)};
        Array<int,AMREX_SPACEDIM> ncells{AMREX_D_DECL(0, 0, 0)};
        GpuArray<int,AMREX_SPACEDIM> periodic{AMREX_D_DECL(0, 0, 0)};
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
        std::string prefix_diagnosticfoldername = "";
        int num_of_digits_in_filenames=6;

        //Diagnostic parameters
        int print_diagnostics=0;
        int ifrigidnodespresent=0;
        bool is_standard_test = true;
        int test_number = 0;	//0->is_standard_test=false,1-->axial vibration of bar, 2-->Dam break, 3--> Elastic collision of disks, 4--> Lateral deflection of beam
        int no_of_rigidbodies_present=0;
        int total_number_of_material_points=0;          //Includes only the number of normal material points(not the rigid particles)
        amrex::Real total_mass_phase_0=0.0;
        amrex::Real total_vol_phase_0=0.0;
        int print_max_min_avg_stress=1;
        int print_max_min_avg_velocity=1;
        int print_total_energy=1;
        int diagnostic_iteration_frequency=10;


        //Rigid body nodes
        Rigid_Bodies *Rb;

        //Parameters for Axial bar test
        amrex::Real axial_bar_E = 100.0;
        amrex::Real axial_bar_rho = 1.0;
        amrex::Real axial_bar_L = 25.0;
        int axial_bar_modenumber = 1;
        amrex::Real axial_bar_v0 = 0.1;

        //Parameters for dam break
        amrex::Real dam_break_H1 = 0.2;
        amrex::Real dam_break_g = 9.81;

        //Parameters for lateral vibration of bar (TVB). The test case is performed for first mode
        amrex::Real tvb_E = 50.0;
        amrex::Real tvb_L = 0.1;
        amrex::Real tvb_rho = 25.0;
        amrex::Real tvb_v0 = 0.0;

        amrex::Real p2g_L = 1.0;
        amrex::Real p2g_f = 1.0;
        int p2g_ncell=0;

        //Parameters for membrane compaction
        amrex::Real mem_compaction_area=1.0;
        amrex::Real mem_compaction_dampcoeff=1.0;
        amrex::Real mem_compaction_L0=1.0;
        int mem_compaction_restoring_force_calc_method = 1;

        //Parameters for spring alone simulation
        amrex::Real spring_alone_area = 0.1;
        amrex::Real spring_alone_length = 0.0;
        amrex::Real spring_alone_E = 0.0;
        amrex::Real spring_alone_exact_deflection = 0.0;
        amrex::Real spring_alone_exact_delta = 0.0;


        //1-->tent function, 2--> Quadratic Spline Shape Function, 3--> Cubic Spline Shape Function
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
        Real total_rigid_mass = 0.0;
        Real total_vol = 0.0;
        Real mem_compaction_vnew=0.0;
        Real mem_compaction_vold=0.0;
        Real mem_compaction_volume=0.0;
        int stress_update_scheme=1;
        int calculate_strain_based_on_delta=0;
        

        //Rigid body params
        int num_of_rigid_bodies=0;

        Vector<int> bclo;
        Vector<int> bchi;

        GpuArray<int,AMREX_SPACEDIM> order_scheme_directional{AMREX_D_DECL(1,1,1)};


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
            Array<int,AMREX_SPACEDIM> periodicarr{AMREX_D_DECL(0, 0, 0)};
            pp.get("is_it_periodic",periodicarr);
            periodic={AMREX_D_DECL(periodicarr[XDIR],periodicarr[YDIR],periodicarr[ZDIR])};
            pp.get("max_grid_size", max_grid_size);



            pp.get("max_steps", maxsteps);
            pp.get("final_time",final_time);

            pp.query("gravity", gravity);
            pp.query("calculate_strain_based_on_delta",calculate_strain_based_on_delta);

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
            pp.query("prefix_diagnosticfoldername",prefix_diagnosticfoldername);

            //Reading diagnostic parameters
            pp.query("is_standard_test",is_standard_test);
            pp.query("test_number",test_number);
            pp.query("axial_bar_E",axial_bar_E);
            pp.query("axial_bar_rho",axial_bar_rho);
            pp.query("axial_bar_L",axial_bar_L);
            pp.query("axial_bar_modenumber",axial_bar_modenumber);
            pp.query("axial_bar_v0",axial_bar_v0);

            //Reading TVB test parameters
            pp.query("tvb_v0",tvb_v0);
            pp.query("tvb_L",tvb_L);
            pp.query("tvb_E",tvb_E);
            pp.query("tvb_rho",tvb_rho);

            pp.query("p2g_L",p2g_L);
            pp.query("p2g_f",p2g_f);
            pp.query("p2g_ncell",p2g_ncell);

            pp.query("mem_compaction_area",mem_compaction_area);
            pp.query("mem_compaction_volume",mem_compaction_volume);
            pp.query("mem_compaction_dampcoeff",mem_compaction_dampcoeff);
            pp.query("mem_compaction_restoring_force_calc_method",mem_compaction_restoring_force_calc_method);

            pp.query("spring_alone_area",spring_alone_area);
            pp.query("spring_alone_length",spring_alone_length);
            pp.query("spring_alone_E",spring_alone_E);

            pp.query("print_max_min_avg_stress",print_max_min_avg_stress);
            pp.query("print_max_min_avg_velocity",print_max_min_avg_velocity);
            pp.query("print_total_energy",print_total_energy);
            pp.query("diagnostic_iteration_frequency",diagnostic_iteration_frequency);


            //Rigid body params
            pp.get("num_of_rigid_bodies", num_of_rigid_bodies);


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

        std::string Convert_arr_string(Array<Real,AMREX_SPACEDIM>arr)
        {
        	std::string str;
        	str="[";
        	for (int i=0;i<AMREX_SPACEDIM-1;i++)
        	{
        		str=str+std::to_string(arr[i])+",";
        	}
        	str=str+std::to_string(arr[AMREX_SPACEDIM-1]);
        	return(str);
        }

        void PrintSimulationParams()
        {
        	//To be completed, Sreejith
        	std::string vect_str;
        	vect_str=Convert_arr_string(plo);
        	amrex::Print()<<"\n\tBackground grid physical extents: "<<vect_str;
        }
        void ThreeDOF_Solver(amrex::Real dt)
        {
        	for (int i=0;i<no_of_rigidbodies_present;i++)
        	{
        		if(Rb[i].position_update_method==1)
        		{
        			for(int j=0;j<AMREX_SPACEDIM;j++)
        			{
        				Rb[i].velocity[j]+=(Rb[i].force_internal[j]+Rb[i].force_external[j])/Rb[i].total_mass*dt;
        			}
        		}
        		else
        		{
        			for(int j=0;j<AMREX_SPACEDIM;j++)
        			{
        				Rb[i].velocity[j]=Rb[i].imposed_velocity[j];
        			}
        		}
        	}
        }

};

//void MPMspecs::PrintSimulationParams()
//{
	//amrex::Print()<<"\n\tBackground grid physical extents: ["<<plo[0]<<","<<"]";
//}
#endif

