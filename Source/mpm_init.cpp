#include <mpm_particle_container.H>
#include <constants.H>

void MPMParticleContainer::InitParticles (const std::string& filename)
{

    // only read the file on the IO proc
    if (ParallelDescriptor::IOProcessor())  
    {
        std::ifstream ifs;
        ifs.open(filename.c_str(), std::ios::in);

        if (!ifs.good())
        {
            amrex::FileOpenFailed(filename);
        }

        int np = -1;
        ifs >> np >> std::ws;

        if ( np == -1 )
        {
            Abort("\nCannot read number of particles from particle file\n");
        }

        const int lev  = 0;
        const int grid = 0;
        const int tile = 0;

        auto& particle_tile = DefineAndReturnParticleTile(lev,grid,tile);
        ParticleType p;

        for (int i = 0; i < np; i++) 
        {
            // Set id and cpu for this particle
            p.id()  = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            // Read from input file
            ifs >> p.idata(intData::phase);
            ifs >> p.pos(0);
            ifs >> p.pos(1);
            ifs >> p.pos(2);
            ifs >> p.rdata(realData::radius);
            ifs >> p.rdata(realData::density);
            ifs >> p.rdata(realData::xvel);
            ifs >> p.rdata(realData::yvel);
            ifs >> p.rdata(realData::zvel);

            // Set other particle properties
            p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three);
            p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
            //amrex::Print()<<"\n Mass = "<<p.rdata(realData::mass);
            p.rdata(realData::jacobian)	   = 1.0;
            //p.rdata(realData::vol_init)	   = 0.0;

            for(int comp=0;comp<NCOMP_TENSOR;comp++)
            {
                p.rdata(realData::strainrate+comp) = zero;
                p.rdata(realData::strain+comp)     = zero;
                p.rdata(realData::stress+comp)     = zero;
            }

            // Add everything to the data structure
            particle_tile.push_back(p);

            if (!ifs.good())
            {
                amrex::Abort("Error initializing particles from Ascii file. \n");
            }
        }
    }
    Redistribute();
}
