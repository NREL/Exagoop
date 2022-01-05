#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <mpm_check_pair.H>
#include <mpm_particle_container.H>
#include <AMReX_PlotFileUtil.H>

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
        mpm_pc.deposit_onto_grid(nodaldata);

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
        WriteSingleLevelPlotfile(pltfile, nodaldata, nodaldata_names, geom, time, 0);
       
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
            {
                mpm_pc.updateNeighbors();
            }
            
            mpm_pc.deposit_onto_grid(nodaldata);
            BL_PROFILE_VAR("MOVE_PART",movepart);
            mpm_pc.moveParticles(dt,specs.gravity);
            BL_PROFILE_VAR_STOP(movepart);

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
                WriteSingleLevelPlotfile(pltfile, nodaldata, nodaldata_names, geom, time, 0);
                
                output_time=zero;
                BL_PROFILE_VAR_STOP(outputs);
            }
            steps++;
        }

        mpm_pc.Redistribute();
        mpm_pc.writeParticles(output_it+1);
        pltfile = amrex::Concatenate("nplt", output_it+1, 5);
        WriteSingleLevelPlotfile(pltfile, nodaldata, nodaldata_names, geom, time, 0);
    }

    amrex::Finalize();
}
