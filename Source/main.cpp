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
        MPMspecs specs;
        specs.read_mpm_specs();

        RealBox real_box;
        for (int n = 0; n < AMREX_SPACEDIM; n++)
        {
            real_box.setLo(n, specs.plo[n]);
            real_box.setHi(n, specs.phi[n]);
        }

        IntVect domain_lo(AMREX_D_DECL(0,0,0));
        IntVect domain_hi(AMREX_D_DECL(specs.ncells[XDIR]-1,
                    specs.ncells[YDIR]-1,
                    specs.ncells[ZDIR]-1));

        const Box domain(domain_lo, domain_hi);

        int coord = 0; //cartesian
        Geometry geom(domain, &real_box, coord, specs.periodic.data());

        BoxArray ba(domain);
        ba.maxSize(specs.max_grid_size);
        DistributionMapping dm(ba);

        const int ng_cells = one;
        MPMParticleContainer mpm_pc(geom, dm, ba, ng_cells);
        mpm_pc.InitParticles(specs.particlefilename);
        
        const BoxArray& nodeba = amrex::convert(ba, IntVect{1,1,1});
        MultiFab nodaldata(nodeba, dm, NUM_STATES, 0);
        nodaldata.setVal(0.0);
        mpm_pc.deposit_onto_grid(nodaldata,specs.gravity,1,0);
        mpm_pc.interpolate_from_grid(nodaldata,0,1);

        int steps=0;
        Real time=zero;
        Real dt=specs.timestep;
        Real output_time=zero;
        Real output_timePrint=zero;
        int output_it=0;

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
       
        amrex::Print() << "Num particles after init is " << mpm_pc.TotalNumberOfParticles() << "\n";


        while((steps < specs.maxsteps) and (time < specs.final_time))
        {
            time += dt;
            output_time += dt;
            output_timePrint += dt;

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
            
            nodaldata.setVal(zero);
            //find mass/vel at nodes
            //update_massvel=1, update_forces=0
            mpm_pc.deposit_onto_grid(nodaldata,specs.gravity,1,0);
            //amrex::Print()<<"Done with mass/vel deposition\n";

            //find strainrate at material points
            //update_vel=0,update_strainrate=1
            mpm_pc.interpolate_from_grid(nodaldata,0,1);
            //amrex::Print()<<"Done with strain-rate interpolation\n";

            //update stress at material points
            mpm_pc.apply_constitutive_model(dt,specs.Youngs_modulus,specs.Poissons_ratio); 
            //amrex::Print()<<"Done with constitutive model\n";

            //update forces at nodes
            //update_massvel=0, update_forces=1
            mpm_pc.deposit_onto_grid(nodaldata,specs.gravity,0,1);
            //amrex::Print()<<"Done with force update\n";

            //update velocity on nodes
            nodal_update(nodaldata,dt);
            //amrex::Print()<<"Done with nodal update\n";

            //impose bcs at nodes
            nodal_bcs(geom,nodaldata,dt);
            //amrex::Print()<<"Done with nodal bcs\n";

            //find velocity at material points
            //update_vel=1,update_strainrate=0
            mpm_pc.interpolate_from_grid(nodaldata,1,0);
            //amrex::Print()<<"Done with velocity i->mp\n";

            //move material points
            mpm_pc.moveParticles(dt);
            //amrex::Print()<<"Done with moving material points\n";
            //amrex::Print()<<"=================================\n";


            if (output_timePrint > specs.screen_output_time)
            {
                Print()<<"step:"<<steps<<"\t"<<"time:"<<time<<"\n";
                output_timePrint=zero;
            }

            if (output_time > specs.write_output_time) 
            {
                BL_PROFILE_VAR("OUTPUT_TIME",outputs);
                Print()<<"writing outputs at step,time:"<<steps<<"\t"<<time<<"\n";
                mpm_pc.Redistribute();
                mpm_pc.fillNeighbors();
                mpm_pc.buildNeighborList(CheckPair());
               
                output_it++;
                mpm_pc.writeParticles(output_it);
                
                pltfile = amrex::Concatenate("nplt", output_it, 5);
                write_plot_file(pltfile,nodaldata,nodaldata_names,geom,ba,dm,time);
                
                output_time=zero;
                BL_PROFILE_VAR_STOP(outputs);
            }
            steps++;
        }

        mpm_pc.Redistribute();
        mpm_pc.writeParticles(output_it+1);
        pltfile = amrex::Concatenate("nplt", output_it+1, 5);
        write_plot_file(pltfile,nodaldata,nodaldata_names,geom,ba,dm,time);
    }

    amrex::Finalize();
}
