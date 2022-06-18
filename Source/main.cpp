#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <mpm_check_pair.H>
#include <mpm_particle_container.H>
#include <AMReX_PlotFileUtil.H>
#include <nodal_update.H>


using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
    	//Reading input files for the simulation
        MPMspecs specs;
        specs.read_mpm_specs();	//Read input file

        //
        int coord = 0; //cartesian
        RealBox real_box;
        for (int n = 0; n < AMREX_SPACEDIM; n++)
        {
            real_box.setLo(n, specs.plo[n]);
            real_box.setHi(n, specs.phi[n]);
        }

        IntVect domain_lo(AMREX_D_DECL(0,0,0));
        IntVect domain_hi(AMREX_D_DECL(specs.ncells[XDIR]-1,specs.ncells[YDIR]-1,specs.ncells[ZDIR]-1));
        const Box domain(domain_lo, domain_hi);

        Geometry geom(domain, &real_box, coord, specs.periodic.data());

        BoxArray ba(domain);
        ba.maxSize(specs.max_grid_size);
        DistributionMapping dm(ba);

        int ng_cells = 1;
        if(specs.order_scheme==3)
        {
        	ng_cells = 3;
        }

        //Initialise particle properties
        MPMParticleContainer mpm_pc(geom, dm, ba, ng_cells);
        mpm_pc.InitParticles(specs.particlefilename,&specs.total_mass,&specs.total_vol);
        
        //Set grid properties
        const BoxArray& nodeba = amrex::convert(ba, IntVect{1,1,1});

        int ng_cells_nodaldata=1;
        if(specs.order_scheme==1)
        {
        	ng_cells_nodaldata=1;
        }
        else if(specs.order_scheme==3)
        {
        	ng_cells_nodaldata=3;
        }

        MultiFab nodaldata(nodeba, dm, NUM_STATES, ng_cells_nodaldata);
        nodaldata.setVal(0.0,ng_cells_nodaldata);

        MultiFab dens_field_data;
        BoxArray dens_ba=ba;
        Box dom_dens = geom.Domain();
        dom_dens.refine(specs.dens_field_gridratio);
        Geometry geom_dens(dom_dens);
        if(specs.dens_field_output)
        {
           dens_ba.refine(specs.dens_field_gridratio);
           const BoxArray& nodal_dens_ba=amrex::convert(dens_ba,IntVect{1,1,1});
           dens_field_data.define(nodal_dens_ba,dm,1,0);
           dens_field_data.setVal(0.0);
        }

        //mpm_pc.fillNeighbors();
        mpm_pc.RedistributeLocal();
                        mpm_pc.fillNeighbors();
                        //mpm_pc.buildNeighborList(CheckPair());
        Real dt=specs.timestep;
        mpm_pc.deposit_onto_grid(nodaldata,specs.gravity,
                                 specs.external_loads_present,
                                 specs.force_slab_lo,
                                 specs.force_slab_hi,
                                 specs.extforce,1,0,specs.mass_tolerance,specs.order_scheme);	//Deposit mass and velocity on node

        mpm_pc.interpolate_mass_from_grid(nodaldata,1);						//Calculate volume of each mp

        mpm_pc.interpolate_from_grid(nodaldata,0,1,specs.order_scheme,specs.alpha_pic_flip);	//Calculate strainrate at each mp
        dt = mpm_pc.Calculate_time_step();
        dt=specs.CFL*dt;
        dt=min(dt,specs.dtmin);


        mpm_pc.apply_constitutive_model(dt,specs.applied_strainrate);

        if(specs.dens_field_output)
        {
           mpm_pc.update_density_field(dens_field_data,specs.dens_field_gridratio,specs.smoothfactor);
        }

        //Quantities for elastic disk collisions
        Real TKE=0.0;
        Real TSE=0.0;
        Real TE=TKE+TSE;

        Real Vmnum=0.0;
        Real Vmex=0.0;
        std::ofstream OutFile;

        /*
        if(amrex::ParallelDescriptor::IOProcessor())
        {

        	std::string FullPathFile = "Energy.out";
        	OutFile.open(FullPathFile.c_str(), std::ios::out);
        	OutFile <<time<<"\t"<<TKE<<"\t"<<TSE<<"\t"<<TE;
        }*/


        //Elastic  collision specific

        int steps=0;
        Real time=zero;
        Real output_time=zero;
        Real output_timePrint=zero;
        int output_it=0;

        if(amrex::ParallelDescriptor::IOProcessor())
                {
                	std::string FullPathFile = "AxialBar.out";
                	OutFile.open(FullPathFile.c_str(), std::ios::out);
                	mpm_pc.CalculateVelocity(Vmnum);
                	Real n = 1;
                	Real beta_n = (2*n-1.0)/2*3.141592/25.0;
                	Real w_n = 10*beta_n;
                	Vmex = 0.1/(beta_n*25.0)*cos(w_n*time);

                	OutFile <<time<<"\t"<<Vmex<<"\t"<<Vmnum;
                }

        mpm_pc.writeParticles(steps);
        amrex::Vector<std::string> nodaldata_names;
        nodaldata_names.push_back("mass");
        nodaldata_names.push_back("vel_x");
        nodaldata_names.push_back("vel_y");
        nodaldata_names.push_back("vel_z");
        nodaldata_names.push_back("force_x");
        nodaldata_names.push_back("force_y");
        nodaldata_names.push_back("force_z");

        std::string pltfile;
        pltfile = amrex::Concatenate("nplt", steps, 5);
        write_plot_file(pltfile,nodaldata,nodaldata_names,geom,ba,dm,time);

        if(specs.dens_field_output)
        {
            pltfile = amrex::Concatenate("dplt", steps, 5);
            write_plot_file(pltfile,dens_field_data,{"density"},geom_dens,dens_ba,dm,time);
        }


        while((steps < specs.maxsteps) and (time < specs.final_time))
        {
        	dt = mpm_pc.Calculate_time_step();
        	dt=specs.CFL*dt;
        	dt=min(dt,specs.dtmin);

            time += dt;
            output_time += dt;
            output_timePrint += dt;

            if (output_timePrint > specs.screen_output_time)
            {
            	Print()<<"step:"<<steps<<"\t"<<"time:"<<time<<" dt = "<<dt<<"\n";
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
            //find mass/vel at nodes
            //update_massvel=1, update_forces=0
            mpm_pc.deposit_onto_grid(nodaldata,specs.gravity,
                                 specs.external_loads_present,
                                 specs.force_slab_lo,
                                 specs.force_slab_hi,
                                 specs.extforce,1,0,specs.mass_tolerance,specs.order_scheme); 		//Update mass and velocity only

            backup_current_velocity(nodaldata);									//Store velocity at time level t to calculate Delta_vel later for flip update
            mpm_pc.deposit_onto_grid(nodaldata,specs.gravity,					// Calculate forces on nodes
                                             specs.external_loads_present,
                                             specs.force_slab_lo,
                                             specs.force_slab_hi,
                                             specs.extforce,0,1,specs.mass_tolerance,specs.order_scheme);
            //update velocity on nodes
            nodal_update(nodaldata,dt,specs.mass_tolerance);

            //impose bcs at nodes
            nodal_bcs(geom,nodaldata,dt);

            //Calculate velocity diff
            store_delta_velocity(nodaldata);

            //Update particle velocity at time t+dt
            mpm_pc.interpolate_from_grid(nodaldata,1,0,specs.order_scheme,specs.alpha_pic_flip);
            mpm_pc.updateNeighbors();

            //Update particle position at t+dt
            mpm_pc.moveParticles(dt);

            if(specs.stress_update_scheme==1)										//MUSL scheme
            {
            	mpm_pc.deposit_onto_grid(nodaldata,specs.gravity,					// Calculate velocity on nodes
            	                                             specs.external_loads_present,
            	                                             specs.force_slab_lo,
            	                                             specs.force_slab_hi,
            	                                             specs.extforce,1,0,specs.mass_tolerance,specs.order_scheme);
            	nodal_bcs(geom,nodaldata,dt);

            }

            //find strainrate at material points at time t+dt
            mpm_pc.interpolate_from_grid(nodaldata,0,1,specs.order_scheme,specs.alpha_pic_flip);
            mpm_pc.updateNeighbors();

            //mpm_pc.move_particles_from_nodevel(nodaldata,dt,1);
            mpm_pc.updatevolume(dt);

            //update stress at material pointsat time t+dt
            if(time<specs.applied_strainrate_time)
            {

                mpm_pc.apply_constitutive_model(dt,
                                                specs.applied_strainrate
												);
            }
            else
            {
                mpm_pc.apply_constitutive_model(dt,
												0.0
												);
            }

            if(specs.dens_field_output)
            {
                mpm_pc.update_density_field(dens_field_data,specs.dens_field_gridratio,specs.smoothfactor);
            }

            mpm_pc.CalculateVelocity(Vmnum);
                            	Real n = 1;
                            	Real beta_n = (2*n-1.0)/2*3.141592/25.0;
                            	Real w_n = 10*beta_n;
                            	Vmex = 0.1/(beta_n*25.0)*cos(w_n*time);

                            	OutFile <<"\n"<<time<<"\t"<<Vmex<<"\t"<<Vmnum;

            if (output_time > specs.write_output_time) 
            {
                BL_PROFILE_VAR("OUTPUT_TIME",outputs);
                Print()<<"writing outputs at step,time:"<<steps<<"\t"<<time<<"\n";
                mpm_pc.Redistribute();
                mpm_pc.fillNeighbors();
                //mpm_pc.buildNeighborList(CheckPair());

                output_it++;
                mpm_pc.writeParticles(output_it);

                pltfile = amrex::Concatenate("nplt", output_it, 5);
                write_plot_file(pltfile,nodaldata,nodaldata_names,geom,ba,dm,time);
            
                if(specs.dens_field_output)
                {
                    pltfile = amrex::Concatenate("dplt", output_it, 5);
                    write_plot_file(pltfile,dens_field_data,{"density"},geom_dens,dens_ba,
                                    dm,time);
                }

                output_time=zero;
                BL_PROFILE_VAR_STOP(outputs);
            }
            steps++;
        }

        mpm_pc.Redistribute();
        mpm_pc.writeParticles(output_it+1);
        
        pltfile = amrex::Concatenate("nplt", output_it+1, 5);
        write_plot_file(pltfile,nodaldata,nodaldata_names,geom,ba,dm,time);
        
        if(specs.dens_field_output)
        {
            pltfile = amrex::Concatenate("dplt", output_it+1, 5);
            write_plot_file(pltfile,dens_field_data,{"density"},geom_dens,dens_ba,dm,time);
        }
    }

    amrex::Finalize();
}
