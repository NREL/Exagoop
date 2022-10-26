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
                                 specs.total_mass,specs.total_vol);
        }
        else
        {   //yli edit for gbhypo
            mpm_pc.InitParticles(specs.autogen_mincoords.data(),specs.autogen_maxcoords.data(),
                                 specs.autogen_vel.data(),specs.autogen_dens,specs.autogen_constmodel,
                                 specs.autogen_E,specs.autogen_nu,
                                 specs.autogen_bulkmod,specs.autogen_Gama_pres,specs.autogen_visc,
                                 specs.autogen_multi_part_per_cell,specs.total_mass,specs.total_vol,specs.initial_void_ratio);
            //end yli edit
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

        //yli edit for gbhypo
        mpm_pc.apply_constitutive_model(specs,dt,specs.applied_strainrate);
        //end yli edit

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
            if(specs.ifrigidnodespresent==1)
            {
            	mpm_pc.deposit_onto_grid_rigidnodesonly(	nodaldata,
            												specs.gravity,
															specs.external_loads_present,
															specs.force_slab_lo,
															specs.force_slab_hi,
															specs.extforce,0,1,
															specs.mass_tolerance,
															specs.order_scheme_directional,
															specs.periodic);
            }

            //update velocity on nodes
            nodal_update(nodaldata,dt,specs.mass_tolerance);

            if(specs.ifrigidnodespresent==1)
            {
            	//The following code statements are not generic. I have coded them just to check the membrane compaction framework (Sreejith, 21 Sept 2022)
            	//Calculate restoring force and velocity
            	amrex::Real vel_piston_new=mpm_pc.GetVelPiston(dt,vel_piston_old);
            	amrex::Print()<<"\n Piston vel = "<<vel_piston_new;

            	vel_piston_old = vel_piston_new;
            	//Detect contact from rigid nodes
            	amrex::Real contact_alpha=specs.mass_tolerance;
            	nodal_detect_contact(nodaldata,contact_alpha,vel_piston_new);
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
                    //yli edit for gbhypo1
            		mpm_pc.apply_constitutive_model(specs,dt,specs.applied_strainrate);
                    //end yli edit
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
                    //yli edit for gbhypo
            		mpm_pc.apply_constitutive_model(specs,dt,0.0);
                    //end yli edit
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
