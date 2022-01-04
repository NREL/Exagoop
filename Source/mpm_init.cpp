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

            p.rdata(realData::fx) = zero;
            p.rdata(realData::fy) = zero;
            p.rdata(realData::fz) = zero;

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

void MPMParticleContainer::addnodalparticles ()
{
    const int lev = 0;
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box = mfi.tilebox();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) 
        {
            amrex::Real  x = plo[XDIR] + iv[XDIR] *dx[0];
            amrex::Real  y = plo[YDIR] + iv[YDIR] *dx[1];
            amrex::Real  z = plo[ZDIR] + iv[ZDIR] *dx[2];

            ParticleType p;
            p.id()  = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();                

            p.pos(XDIR) = x;
            p.pos(YDIR) = y;
            p.pos(ZDIR) = z;

            p.idata(intData::phase) = NPOINT;
            p.rdata(realData::radius)  = std::pow(dx[0]*dx[1]*dx[2],0.333333);
            p.rdata(realData::xvel)    = 0.0;
            p.rdata(realData::yvel)    = 0.0;
            p.rdata(realData::zvel)    = 0.0;
            p.rdata(realData::fx)      = 0.0;
            p.rdata(realData::fy)      = 0.0;
            p.rdata(realData::fz)      = 0.0;
            p.rdata(realData::volume)  = fourbythree*PI*pow(p.rdata(realData::radius),three);
            p.rdata(realData::mass)    = 0.0;
            p.rdata(realData::density) = 0.0;

            particle_tile.push_back(p);
        }
    }

    Redistribute();
}
