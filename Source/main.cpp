#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <mpm_check_pair.H>
#include <mpm_particle_container.H>
#include <constants.H>

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

        int steps=0;
        Real time=zero;
        Real dt=specs.timestep;
        Real output_time=zero;
        Real output_timePrint=zero;
        int output_it=0;
        mpm_pc.writeParticles(steps);
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
                output_time=zero;
                BL_PROFILE_VAR_STOP(outputs);
            }
            steps++;
        }

        mpm_pc.Redistribute();
        mpm_pc.writeParticles(output_it+1);
    }

    amrex::Finalize();
}
