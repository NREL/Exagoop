#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <mpm_check_pair.H>
#include <mpm_particle_container.H>
#include <AMReX_PlotFileUtil.H>
#include <nodal_data_ops.H>
#include <mpm_eb.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
    	//Initializing and reading input file for the simulation
        MPMspecs specs;
        specs.read_mpm_specs();

        //Declaring solver variables
        int steps=0;
        Real dt;
        Real time = 0.0;
        int output_it=0;
        std::string pltfile;
        Real output_time=zero;
        Real output_timePrint=zero;

        Array <int,AMREX_SPACEDIM> order_surface_integral={3,3,3};

        //Check if max_grid_size==1. If yes, Abort.
        if(specs.max_grid_size==1)
        {
        	amrex::Abort("\nMax grid size should be greater than or equal to two");
        }

        int coord = 0; //cartesian
        RealBox real_box;
        for (int n = 0; n < AMREX_SPACEDIM; n++)
        {
            real_box.setLo(n, specs.plo[n]);
            real_box.setHi(n, specs.phi[n]);
        }

        IntVect domain_lo(AMREX_D_DECL(0,0,0));
        IntVect domain_hi(AMREX_D_DECL(	specs.ncells[XDIR]-1,
                    					specs.ncells[YDIR]-1,
										specs.ncells[ZDIR]-1));
        const Box domain(domain_lo, domain_hi);
        Geometry geom(domain, &real_box, coord, specs.periodic.data());

        //Create box array and chunking it
        BoxArray ba(domain);
        ba.maxSize(specs.max_grid_size);
        DistributionMapping dm(ba);

        //Defining number of ghost cells for particle data
        int ng_cells = 1;
        if(specs.order_scheme==3)
        {
        	ng_cells = 2;
        }
        mpm_ebtools::init_eb(geom,ba,dm);


        //Initialise particle properties
        MPMParticleContainer mpm_pc(geom, dm, ba, ng_cells);

        if(specs.restart_checkfile !="")	//Restart from checkpoint solution
        {
        	mpm_pc.readCheckpointFile(specs.restart_checkfile, steps,time,output_it);
        	Print()<<"\nRestarting from checkpoint file: "<<specs.restart_checkfile;
        }
        else if(!specs.use_autogen)
        {
            mpm_pc.InitParticles(specs.particlefilename,
                                 specs.total_mass,specs.total_vol,specs.total_rigid_mass);
        }
        else
        {
            mpm_pc.InitParticles(specs.autogen_mincoords.data(),specs.autogen_maxcoords.data(),
                                 specs.autogen_vel.data(),specs.autogen_dens,specs.autogen_constmodel,
                                 specs.autogen_E,specs.autogen_nu,
                                 specs.autogen_bulkmod,specs.autogen_Gama_pres,specs.autogen_visc,
                                 specs.autogen_multi_part_per_cell,specs.total_mass,specs.total_vol);
        }

        specs.ifrigidnodespresent = mpm_pc.checkifrigidnodespresent();

        if(mpm_ebtools::using_levelset_geometry)
        {
            mpm_pc.removeParticlesInsideEB();
        }

        //Set background grid properties
        const BoxArray& nodeba = amrex::convert(ba, IntVect{1,1,1});

        int ng_cells_nodaldata=1;
        if(specs.order_scheme==1)
        {
            ng_cells_nodaldata=1;
        }
        else if(specs.order_scheme==3)
        {
            ng_cells_nodaldata=3;

            specs.order_scheme_directional[XDIR] = ((specs.periodic[XDIR]==0)?((specs.ncells[XDIR]<5)?1:3):((specs.ncells[XDIR]<3)?1:3));
            specs.order_scheme_directional[YDIR] = ((specs.periodic[YDIR]==0)?((specs.ncells[YDIR]<5)?1:3):((specs.ncells[YDIR]<3)?1:3));
            specs.order_scheme_directional[ZDIR] = ((specs.periodic[ZDIR]==0)?((specs.ncells[ZDIR]<5)?1:3):((specs.ncells[ZDIR]<3)?1:3));

            if(specs.order_scheme_directional[XDIR]==1 and specs.order_scheme_directional[YDIR]==1 and specs.order_scheme_directional[ZDIR]==1 )
            {
            	amrex::Print()<<"\nWarning! Number of cells in all directions do not qualify for cubic-spline shape functions. Hence calculating using linear hat shape functions in all directions";
            }

            //Make sure that none of the boxes that use spline function are of size of 1.
            //For example if ncell=5 and max_grid_size = 2,we get boxes of {2,2,1}. I (Sreejith) noticed that when the box size is one
            //ghost particles are not placed correctly.

            for(int box_index=0;box_index<ba.size();box_index++)
            {
            	for(int dim=0;dim<AMREX_SPACEDIM;dim++)
            	{
            		if(ba[box_index].size()[dim]==1 and specs.order_scheme_directional[dim]==3)
            		{
            		   	amrex::Abort("\nError! Box cannot be of size =1. Please change max_grid_size value to make sure all boxes have size>1 when using splines");
            		}
            	}
            }
        }
        else
        {
            amrex::Abort("\nOrder scheme not implemented yet. Please use order_scheme=1 or order_scheme=3 in the input file.\n");
        }

        MultiFab nodaldata(nodeba, dm, NUM_STATES, ng_cells_nodaldata);
        nodaldata.setVal(0.0,ng_cells_nodaldata);

        MultiFab dens_field_data;
        BoxArray dens_ba=ba;
        Box dom_dens = geom.Domain();
        dom_dens.refine(specs.dens_field_gridratio);
        Geometry geom_dens(dom_dens);
        int ng_dens=3;
        if(specs.dens_field_output)
        {
            dens_ba.refine(specs.dens_field_gridratio);
            dens_field_data.define(dens_ba,dm,1,ng_dens);
            dens_field_data.setVal(0.0,ng_dens);
        }

        //mpm_pc.fillNeighbors();
        mpm_pc.RedistributeLocal();
        mpm_pc.fillNeighbors();

        //Calculate time step
        dt 	= (specs.fixed_timestep==1)?specs.timestep:mpm_pc.Calculate_time_step(specs.CFL,specs.dt_max_limit,specs.dt_min_limit);

        //Deposit mass and velocity on node
        amrex::Print()<<"\nOrder = "<<specs.order_scheme_directional;

        mpm_pc.deposit_onto_grid(nodaldata,
        		 	 	 	 	 specs.gravity,
                                 specs.external_loads_present,
                                 specs.force_slab_lo,
                                 specs.force_slab_hi,
                                 specs.extforce,1,0,specs.mass_tolerance,
								 specs.order_scheme_directional,
								 specs.periodic);

        //Calculate strainrate at each material point
       mpm_pc.interpolate_from_grid(nodaldata,
        							0,
									1,
									specs.order_scheme_directional,
									specs.periodic,
									specs.alpha_pic_flip,
									dt);	//Calculate strainrate at each mp


        mpm_pc.apply_constitutive_model(dt,specs.applied_strainrate);

        if(specs.dens_field_output)
        {
            mpm_pc.update_density_field(dens_field_data,specs.dens_field_gridratio,specs.smoothfactor);
        }

        //Quantities for elastic disk collisions

        if(specs.print_diagnostics)
        {
        	Real TKE=0.0;
        	Real TSE=0.0;
        	Real TE=TKE+TSE;

        	if(specs.is_standard_test)
        	{
        		Real Vmnum=0.0;
        		Real Vmex=0.0;
        		Real Xwf;
        		Real err = 0.0;
        		Real Fy_bottom = 0.0;
        		Real Fy_top = 0.0;

        		switch(specs.test_number)
        		{
        		case(1):	//Axial vibration of continuum bar
							mpm_pc.CalculateVelocity(Vmnum);
        					Vmex = mpm_pc.CalculateExactVelocity(specs.axial_bar_modenumber,specs.axial_bar_E,specs.axial_bar_rho,specs.axial_bar_v0,specs.axial_bar_L,time);
        					PrintToFile("AxialBarVel.out")<<time<<"\t"<<Vmex<<"\t"<<Vmnum<<"\n";
        					mpm_pc.CalculateEnergies(TKE,TSE);
        					TE=TKE+TSE;
        					PrintToFile("AxialBarEnergy.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
        					break;
        		case(2):	//Dam break
							mpm_pc.FindWaterFront(Xwf);
        					PrintToFile("DamBreakWaterfront.out")<<time/sqrt(specs.dam_break_H1/specs.dam_break_g)<<"\t"<<Xwf/specs.dam_break_H1<<"\n";
        					break;
        		case(3):	//Elastic collision of disks
							mpm_pc.CalculateEnergies(TKE,TSE);
        					TE=TKE+TSE;
        		        	PrintToFile("ElasticDiskCollisionEnergy.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
        		        	break;
        		case(4):	//Static deflection of beam under gravity
							mpm_pc.CalculateVelocityCantilever(Vmnum);
				            Vmex = 0.0;
				            PrintToFile("CantileverVel.out")<<time<<"\t"<<Vmex<<"\t"<<Vmnum<<"\n";
				            mpm_pc.CalculateEnergies(TKE,TSE);
				            TE=TKE+TSE;
				            PrintToFile("CantileverEnergy.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
        					break;
        		case(5):	//Transverse vibration of a bar
							mpm_pc.CalculateErrorTVB(specs.tvb_E,specs.tvb_v0,specs.tvb_L,specs.tvb_rho,err);
        					PrintToFile("TVB_Error.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
        					break;
        		case(6):    //Check function reconstruction and convergence
							mpm_pc.CalculateErrorP2G(nodaldata,specs.p2g_L,specs.p2g_f,specs.p2g_ncell);
		        			break;
        		case(7):	//Checks the weight of a block of elastic solid. Used to validate functionality to evaluate surface forces
							mpm_pc.deposit_onto_grid(nodaldata,
		        		 	 	 	 	 specs.gravity,
		                                 specs.external_loads_present,
		                                 specs.force_slab_lo,
		                                 specs.force_slab_hi,
		                                 specs.extforce,0,2,specs.mass_tolerance,
										 specs.order_scheme_directional,
										 specs.periodic);
        					CalculateSurfaceIntegralOnBG(geom, nodaldata,STRESS_INDEX,err);
        					PrintToFile("Weight.out")<<time<<"\t"<<err<<"\t"<<-specs.total_mass*9.81<<"\n";
        					break;
        		case(8): 	CalculateInterpolationError(geom, nodaldata,STRESS_INDEX);
        				    break;
        		case(9):    /*mpm_pc.deposit_onto_grid(nodaldata,
													 specs.gravity,
													 specs.external_loads_present,
													 specs.force_slab_lo,
													 specs.force_slab_hi,
													 specs.extforce,
													 0,
													 2,
													 specs.mass_tolerance,
													 order_surface_integral,
													 specs.periodic);
        					CalculateSurfaceIntegralOnBG(geom, nodaldata,STRESS_INDEX,Fy_bottom);
        					Fy_bottom=-1*Fy_bottom;
        					mpm_pc.CalculateSurfaceIntegralTop(specs.gravity,Fy_top,Fy_bottom);
        					specs.mem_compaction_vnew = mpm_pc.GetVelPiston(dt,specs.mem_compaction_vold,specs.gravity,Fy_top);
        					nodal_detect_contact(nodaldata,specs.mass_tolerance,specs.mem_compaction_vnew);
        					specs.mem_compaction_vold=specs.mem_compaction_vnew;*/
        					specs.mem_compaction_L0=specs.total_vol/specs.mem_compaction_area;
        					break;
        		case(10):   //Spring alone deflection problem.
							//Calculate the eaxct steady state deflection
							specs.spring_alone_exact_deflection = specs.spring_alone_length-specs.total_mass*fabs(specs.gravity[YDIR])/(2.0*specs.spring_alone_E*specs.spring_alone_area/specs.spring_alone_length);
        					specs.spring_alone_exact_delta = specs.spring_alone_length-mpm_pc.GetPosSpring();
        					amrex::Print()<<"\n"<<specs.total_mass<<" "<<specs.gravity[YDIR]<<" "<<specs.spring_alone_length<<" "<<specs.spring_alone_area<<" "<<specs.spring_alone_E;
        				    break;
        		default:	//
        					amrex::Abort("\nSorry. The test number does not exist");

        		}
        	}
        	else
        	{
        		mpm_pc.CalculateEnergies(TKE,TSE);
        		TE=TKE+TSE;
        		PrintToFile("Energy.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
        	}
        }

        amrex::Vector<std::string> nodaldata_names;
        nodaldata_names.push_back("mass");
        nodaldata_names.push_back("vel_x");
        nodaldata_names.push_back("vel_y");
        nodaldata_names.push_back("vel_z");
        nodaldata_names.push_back("force_x");
        nodaldata_names.push_back("force_y");
        nodaldata_names.push_back("force_z");
        nodaldata_names.push_back("delta_velx");
        nodaldata_names.push_back("delta_vely");
        nodaldata_names.push_back("delta_velz");
        nodaldata_names.push_back("mass_old");
        nodaldata_names.push_back("VELX_RIGID_INDEX");
        nodaldata_names.push_back("VELY_RIGID_INDEX");
        nodaldata_names.push_back("VELZ_RIGID_INDEX");
        nodaldata_names.push_back("MASS_RIGID_INDEX");
        nodaldata_names.push_back("STRESS_INDEX");

        if(specs.restart_checkfile =="")
        {
        	mpm_pc.writeParticles(specs.prefix_particlefilename, specs.num_of_digits_in_filenames, steps);

        	pltfile = amrex::Concatenate(specs.prefix_gridfilename, steps,specs.num_of_digits_in_filenames);
        	write_plot_file(pltfile,nodaldata,nodaldata_names,geom,ba,dm,time);

        	if(specs.dens_field_output)
        	{
        		pltfile = amrex::Concatenate(specs.prefix_densityfilename, steps, specs.num_of_digits_in_filenames);
        		WriteSingleLevelPlotfile(pltfile, dens_field_data, {"density"}, geom_dens, time, 0);
        	}
        }

        amrex::Print()<<"\nNumber of particles in the simulation:"<<mpm_pc.TotalNumberOfParticles()<<"\n";
        amrex::Real vel_piston_old=0.0;

        while((steps < specs.maxsteps) and (time < specs.final_time))
        {
        	dt 	= (specs.fixed_timestep==1)?specs.timestep:mpm_pc.Calculate_time_step(specs.CFL,specs.dt_max_limit,specs.dt_min_limit);

            time += dt;
            output_time += dt;
            output_timePrint += dt;
            steps++;

            if (output_timePrint >= specs.screen_output_time)
            {
                Print()<<"\nIteration: "<<std::setw(10)<<steps<<",\t"<<"Time: "<<std::fixed<<std::setprecision(10)<<time<<",\tDt = "<<std::scientific<<std::setprecision(5)<<dt;
                output_timePrint=zero;
            }

            if (steps % specs.num_redist == 0)
            {
                mpm_pc.RedistributeLocal();
                mpm_pc.fillNeighbors();
                mpm_pc.buildNeighborList(CheckPair());
            }
            else 
            {
                mpm_pc.updateNeighbors();
            }

            nodaldata.setVal(zero,ng_cells_nodaldata);

            //update_massvel=1, update_forces=0
            mpm_pc.deposit_onto_grid(	nodaldata,
            							specs.gravity,
										specs.external_loads_present,
										specs.force_slab_lo,
										specs.force_slab_hi,
										specs.extforce,
										1,
										0,
										specs.mass_tolerance,
										specs.order_scheme_directional,
										specs.periodic);

            //Store node velocity at time level t to calculate Delta_vel later for flip update
            backup_current_velocity(nodaldata);									

            // Calculate forces on nodes
            mpm_pc.deposit_onto_grid(	nodaldata,
            							specs.gravity,
										specs.external_loads_present,
										specs.force_slab_lo,
										specs.force_slab_hi,
										specs.extforce,0,1,
										specs.mass_tolerance,
										specs.order_scheme_directional,
										specs.periodic);

            //Calculate mass and velocity from rigid nodes
            if(specs.is_standard_test==1 and specs.test_number ==9 and specs.ifrigidnodespresent==1)
            {
            	mpm_pc.deposit_onto_grid_rigidnodesonly(	nodaldata,
            												specs.gravity,
															specs.external_loads_present,
															specs.force_slab_lo,
															specs.force_slab_hi,
															specs.extforce,1,1,
															specs.mass_tolerance,
															specs.order_scheme_directional,
															specs.periodic);
            }

            //update velocity on nodes
            nodal_update(nodaldata,dt,specs.mass_tolerance);

            if(specs.is_standard_test==1 and specs.test_number ==9 and specs.ifrigidnodespresent==1)
            {
            	Real Fy_top=0.0;
            	Real Fy_bottom=0.0;
            	Real ymin=0.0;

            	if(specs.mem_compaction_restoring_force_calc_method==1)
            	{
            		CalculateSurfaceIntegralOnBG(geom, nodaldata,STRESS_INDEX,Fy_bottom);
            		mpm_pc.CalculateSurfaceIntegralTop(specs.gravity,Fy_top,Fy_bottom);
            	}
            	else
            	{
            		Fy_top=mpm_pc.CalculateEffectiveSpringConstant(specs.mem_compaction_area,specs.mem_compaction_L0);

            	}

            	specs.mem_compaction_vnew = mpm_pc.GetVelPiston(dt,specs.mem_compaction_vold,specs.gravity,Fy_top);
            	nodal_detect_contact(nodaldata,specs.mass_tolerance,specs.mem_compaction_vnew);
            	specs.mem_compaction_vold=specs.mem_compaction_vnew;
            	ymin = mpm_pc.GetPosPiston();
            	PrintToFile("Spring.out")<<time<<"\t"<<ymin<<"\n";
            }


            //impose bcs at nodes
            nodal_bcs(	geom,nodaldata,
            			specs.bclo.data(),
						specs.bchi.data(),
						specs.wall_mu_lo.data(),
						specs.wall_mu_hi.data(),
						specs.wall_vel_lo.data(),
						specs.wall_vel_hi.data(),
						dt);

            if(mpm_ebtools::using_levelset_geometry)
            {
                nodal_levelset_bcs(nodaldata,geom,dt,specs.levelset_bc,
                        specs.levelset_wall_mu);
            }

            //Calculate velocity diff
            store_delta_velocity(nodaldata);

            //Update particle velocity at time t+dt
            mpm_pc.updateNeighbors();
            mpm_pc.interpolate_from_grid(	nodaldata,
            								1,
											0,
											specs.order_scheme_directional,
											specs.periodic,
											specs.alpha_pic_flip,
											dt);
            mpm_pc.updateNeighbors();

            //Update particle position at t+dt
            mpm_pc.moveParticles(	dt,
            						specs.bclo.data(),specs.bchi.data(),
									specs.levelset_bc,
									specs.wall_mu_lo.data(),
									specs.wall_mu_hi.data(),
									specs.wall_vel_lo.data(),
									specs.wall_vel_hi.data(),
									specs.levelset_wall_mu);

            if(specs.stress_update_scheme==1)										
            {
                //MUSL scheme
                // Calculate velocity on nodes
                mpm_pc.deposit_onto_grid(	nodaldata,
                							specs.gravity,
											specs.external_loads_present,
											specs.force_slab_lo,
											specs.force_slab_hi,
											specs.extforce,
											1,
											0,
											specs.mass_tolerance,
											specs.order_scheme_directional,
											specs.periodic);

                nodal_bcs(	geom,
                			nodaldata,
							specs.bclo.data(),
							specs.bchi.data(),
							specs.wall_mu_lo.data(),
							specs.wall_mu_hi.data(),
							specs.wall_vel_lo.data(),
							specs.wall_vel_hi.data(),
							dt);
                //nodal_bcs(	geom, nodaldata, dt);
                
                if(mpm_ebtools::using_levelset_geometry)
                {
                    nodal_levelset_bcs(	nodaldata,geom,
                    					dt,
										specs.levelset_bc,
										specs.levelset_wall_mu);
                }
            }

            //find strainrate at material points at time t+dt
            mpm_pc.interpolate_from_grid(nodaldata,0,1,specs.order_scheme_directional,specs.periodic,specs.alpha_pic_flip,dt);
            mpm_pc.updateNeighbors();

            //mpm_pc.move_particles_from_nodevel(nodaldata,dt,specs.bclo.data(),specs.bchi.data(),1);
            mpm_pc.updateVolume(dt);

            //update stress at material pointsat time t+dt
            if(time<specs.applied_strainrate_time)
            {

            	if(specs.calculate_strain_based_on_delta==1)
            	{
            		mpm_pc.apply_constitutive_model_delta(dt,specs.applied_strainrate);
            	}
            	else
            	{
            		mpm_pc.apply_constitutive_model(dt,specs.applied_strainrate);

            	}
            }
            else
            {
            	if(specs.calculate_strain_based_on_delta==1)
            	{
            		mpm_pc.apply_constitutive_model_delta(dt,0.0);
            	}
            	else
            	{
            		mpm_pc.apply_constitutive_model(dt,0.0);
            	}

            }

            if(specs.dens_field_output)
            {
                mpm_pc.update_density_field(	dens_field_data,
                        						specs.dens_field_gridratio,
												specs.smoothfactor);
            }

            if(specs.print_diagnostics)
            {
            	Real TKE=0.0;
            	Real TSE=0.0;
            	Real TE=TKE+TSE;

            	if(specs.is_standard_test)
            	{
            		Real Vmnum=0.0;
            		Real Vmex=0.0;
            		Real Xwf;
            		Real err;
            		Real Fy_top=0.0;
            		Real Fy_bottom=0.0;
            		Real ymin=0.0;
            		Real ymax = 0.0;


            		switch(specs.test_number)
            		{
            			case(1):	//Axial vibration of continuum bar
									mpm_pc.CalculateVelocity(Vmnum);
                    				Vmex = mpm_pc.CalculateExactVelocity(specs.axial_bar_modenumber,specs.axial_bar_E,specs.axial_bar_rho,specs.axial_bar_v0,specs.axial_bar_L,time);
                    				PrintToFile("AxialBarVel.out")<<time<<"\t"<<Vmex<<"\t"<<Vmnum<<"\n";
                    				mpm_pc.CalculateEnergies(TKE,TSE);
                    				TE=TKE+TSE;
                    				PrintToFile("AxialBarEnergy.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
                    				break;
                    	case(2):	//Dam break
            						mpm_pc.FindWaterFront(Xwf);
                    				PrintToFile("DamBreakWaterfront.out")<<time/sqrt(specs.dam_break_H1/specs.dam_break_g)<<"\t"<<Xwf/specs.dam_break_H1<<"\n";
                    				break;
                    	case(3):	//Elastic collision of disks
            						mpm_pc.CalculateEnergies(TKE,TSE);
                    				TE=TKE+TSE;
                    		       	PrintToFile("ElasticDiskCollisionEnergy.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
                    		       	break;
                    	case(4):	//Static deflection of a beam under gravity
									mpm_pc.CalculateVelocityCantilever(Vmnum);
		                    		Vmex = 0.0;
		                    		PrintToFile("CantileverVel.out")<<time<<"\t"<<Vmex<<"\t"<<Vmnum<<"\n";
		                    		mpm_pc.CalculateEnergies(TKE,TSE);
		                    		TE=TKE+TSE;
		                    		PrintToFile("CantileverEnergy.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
                    				break;
                    	case(5):	//Transverse vibration of a bar
                    				mpm_pc.CalculateErrorTVB(specs.tvb_E,specs.tvb_v0,specs.tvb_L,specs.tvb_rho,err);
                    	        	PrintToFile("TVB_Error.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
                    	        	break;
                    	case(6):    //Check function reconstruction and convergence
									mpm_pc.CalculateErrorP2G(nodaldata,specs.p2g_L,specs.p2g_f,specs.p2g_ncell);
                    			    break;
                    	case(7):	//Checks the weight of a block of elastic solid. Used to validate functionality to evaluate surface forces
									mpm_pc.deposit_onto_grid(nodaldata,
				        		 	 	 	 	 specs.gravity,
				                                 specs.external_loads_present,
				                                 specs.force_slab_lo,
				                                 specs.force_slab_hi,
				                                 specs.extforce,0,2,specs.mass_tolerance,
												 order_surface_integral,
												 specs.periodic);
									CalculateSurfaceIntegralOnBG(geom, nodaldata,STRESS_INDEX,err);
									PrintToFile("Weight.out")<<time<<"\t"<<err<<"\t"<<-specs.total_mass*9.81<<"\n";
                    				break;
                    	case(8):    CalculateInterpolationError(geom, nodaldata,STRESS_INDEX);
                    			    break;
                    	case(9):	/*mpm_pc.deposit_onto_grid(nodaldata,
                    										specs.gravity,
															specs.external_loads_present,
															specs.force_slab_lo,
															specs.force_slab_hi,
															specs.extforce,
															0,
															2,
															specs.mass_tolerance,
															order_surface_integral,
															specs.periodic);
                    				CalculateSurfaceIntegralOnBG(geom, nodaldata,STRESS_INDEX,Fy_bottom);
                    				Fy_bottom=-1*Fy_bottom;
                    				amrex::Print()<<"\n Fy = "<<Fy_bottom;
                    				mpm_pc.CalculateSurfaceIntegralTop(specs.gravity,Fy_top,Fy_bottom);
                    				specs.mem_compaction_vnew = mpm_pc.GetVelPiston(dt,specs.mem_compaction_vold,specs.gravity,Fy_top);
                    				nodal_detect_contact(nodaldata,specs.mass_tolerance,specs.mem_compaction_vnew);
                    				specs.mem_compaction_vold=specs.mem_compaction_vnew;
                    				ymin = mpm_pc.GetPosPiston();
                    				PrintToFile("Spring.out")<<time<<"\t"<<ymin<<"\n";*/
                    				break;
                    	case(10):	//Get oscillations of a single spring under self weight
									ymax = mpm_pc.GetPosSpring()+specs.spring_alone_exact_delta;
                    				PrintToFile("Spring.out")<<time<<"\t"<<ymax<<"\t"<<specs.spring_alone_exact_deflection<<"\n";
                    				break;
                    	default:	//
                    				amrex::Print()<<"\nTest number = "<<specs.test_number;
                    				amrex::Abort("\nSorry. The test number does not exist");

            		}
            	}
            	else
                {
            		mpm_pc.CalculateEnergies(TKE,TSE);
                    TE=TKE+TSE;
                    PrintToFile("Energy.out")<<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE<<"\n";
                 }
            }

            if (fabs(output_time-specs.write_output_time)<dt*0.5)
            {
                BL_PROFILE_VAR("OUTPUT_TIME",outputs);
                Print()<<"\nWriting outputs at step,time:"<<steps<<"\t"<<time;
                mpm_pc.Redistribute();
                mpm_pc.fillNeighbors();

                output_it++;
                mpm_pc.writeParticles(specs.prefix_particlefilename, specs.num_of_digits_in_filenames, output_it);

                pltfile = amrex::Concatenate(specs.prefix_gridfilename, output_it,specs.num_of_digits_in_filenames );
                write_plot_file(pltfile,nodaldata,nodaldata_names,geom,ba,dm,time);

                if(specs.test_number==5)
                {
                	mpm_pc.WriteDeflectionTVB(specs.tvb_E,specs.tvb_v0,specs.tvb_L,specs.tvb_rho,time,output_it);
                }

                if(specs.dens_field_output)
                {
                    pltfile = amrex::Concatenate(specs.prefix_densityfilename, output_it, specs.num_of_digits_in_filenames);
                    WriteSingleLevelPlotfile(pltfile, dens_field_data, {"density"}, geom_dens, time, 0);
                }

                output_time=zero;
                BL_PROFILE_VAR_STOP(outputs);

                mpm_pc.writeCheckpointFile(specs.prefix_checkpointfilename, specs.num_of_digits_in_filenames, time,steps,output_it);
            }

        }

        mpm_pc.Redistribute();
        //mpm_pc.fillNeighbors();
        mpm_pc.writeParticles(specs.prefix_particlefilename, specs.num_of_digits_in_filenames,output_it+1);

        pltfile = amrex::Concatenate(specs.prefix_gridfilename, output_it+1, specs.num_of_digits_in_filenames);
        write_plot_file(pltfile,nodaldata,nodaldata_names,geom,ba,dm,time);
        if(specs.print_diagnostics and specs.is_standard_test and specs.test_number==4)
        {
        	mpm_pc.WriteDeflectionCantilever();
        }

        if(specs.dens_field_output)
        {
            pltfile = amrex::Concatenate(specs.prefix_densityfilename, output_it+1, specs.num_of_digits_in_filenames);
            WriteSingleLevelPlotfile(pltfile, dens_field_data, {"density"}, geom_dens, time, 0);
        }
    }

    amrex::Finalize();
}
