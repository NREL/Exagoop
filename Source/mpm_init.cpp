#include <mpm_particle_container.H>
#include <constants.H>
#include <mpm_eb.H>

void MPMParticleContainer::InitParticles (const std::string& filename,
                                          Real &total_mass,Real &total_vol)
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

        total_mass=0.0;
        total_vol=0.0;

        auto& particle_tile = DefineAndReturnParticleTile(lev,grid,tile);
        Gpu::HostVector<ParticleType> host_particles;

        for (int i = 0; i < np; i++) 
        {
            ParticleType p;
            int ph;
       	    amrex::Real junk;
            // Set id and cpu for this particle
            p.id()  = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            // Read from input file
            //ifs >> junk;
            ifs >> p.idata(intData::phase);
            p.idata(intData::phase) = 0;
            ifs >> p.pos(0);
            ifs >> p.pos(1);
            ifs >> p.pos(2);
            ifs >> p.rdata(realData::radius);
            ifs >> p.rdata(realData::density);
            ifs >> p.rdata(realData::xvel);
            ifs >> p.rdata(realData::yvel);
            ifs >> p.rdata(realData::zvel);
            ifs >> p.idata(intData::constitutive_model);		

            if(p.idata(intData::constitutive_model)==0)	//Elastic solid
            {
            	ifs >> p.rdata(realData::E);
            	ifs >> p.rdata(realData::nu);
            	p.rdata(realData::Bulk_modulus)=0.0;
            	p.rdata(realData::Gama_pressure)=0.0;
            	p.rdata(realData::Dynamic_viscosity)=0.0;
            }
            else if(p.idata(intData::constitutive_model)==1)
            {
            	p.rdata(realData::E)=0.0;
            	p.rdata(realData::nu)=0.0;
            	ifs >> p.rdata(realData::Bulk_modulus);
            	ifs >> p.rdata(realData::Gama_pressure);
            	ifs >> p.rdata(realData::Dynamic_viscosity);
            }
            else
            {
            	amrex::Abort("\nIncorrect constitutive model. Please check your particle file");
            }

            // Set other particle properties
            p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three);	
            //This is the right mass of each particle. Make sure the radius is entered correctly while generating the particle file
            p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);

            total_mass +=p.rdata(realData::mass);
            total_vol +=p.rdata(realData::volume);
            p.rdata(realData::jacobian)	   = 1.0;
            p.rdata(realData::vol_init)	   = p.rdata(realData::volume);
            p.rdata(realData::pressure)    = 0.0;

            for(int comp=0;comp<NCOMP_FULLTENSOR;comp++)
            {
            	p.rdata(realData::deformation_gradient+comp) = 0.0;
            }
            p.rdata(realData::deformation_gradient+0) = 1.0;
            p.rdata(realData::deformation_gradient+4) = 1.0;
            p.rdata(realData::deformation_gradient+8) = 1.0;

            for(int comp=0;comp<NCOMP_TENSOR;comp++)
            {
                p.rdata(realData::strainrate+comp) = zero;
                p.rdata(realData::strain+comp)     = zero;
                p.rdata(realData::stress+comp)     = zero;
            }
            
            host_particles.push_back(p);

            if (!ifs.good())
            {
                amrex::Abort("Error initializing particles from Ascii file. \n");
            }
        }
        
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }
    Redistribute();
}

void MPMParticleContainer::InitParticles (Real mincoords[AMREX_SPACEDIM],Real maxcoords[AMREX_SPACEDIM], 
        Real vel[AMREX_SPACEDIM],
        Real dens, int constmodel, 
        Real E, Real nu,Real bulkmod, Real Gama_pres,Real visc,
        int do_multi_part_per_cell,Real &total_mass,Real &total_vol)
{
    int lev = 0;
    Real x,y,z,x0,y0,z0;

    Real dx = Geom(lev).CellSize(0);
    Real dy = Geom(lev).CellSize(1);
    Real dz = Geom(lev).CellSize(2);
    const Real* plo = Geom(lev).ProbLo();

    total_mass=0.0;
    total_vol=0.0;

    //std::mt19937 mt(0451);
    //std::uniform_real_distribution<double> dist(0.4, 0.6);

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) 
    {

        const Box& tile_box = mfi.tilebox();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
        
        Gpu::HostVector<ParticleType> host_particles;

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) 
        {
            if(do_multi_part_per_cell == 0)
            {
                x = plo[XDIR] + (iv[XDIR] + half)*dx;
                y = plo[YDIR] + (iv[YDIR] + half)*dy;
                z = plo[ZDIR] + (iv[ZDIR] + half)*dz;

                if(x>=mincoords[XDIR] && x<=maxcoords[XDIR] &&
                   y>=mincoords[YDIR] && y<=maxcoords[YDIR] &&
                   z>=mincoords[ZDIR] && z<=maxcoords[ZDIR])
                {
                    ParticleType p = generate_particle(x,y,z,vel,
                            dens,dx*dy*dz,constmodel,
                            E,nu,bulkmod,Gama_pres,visc);

                    total_mass += p.rdata(realData::mass);
                    total_vol += p.rdata(realData::volume);

                    host_particles.push_back(p);
                }
            }
            else
            {
                x0 = plo[XDIR]+iv[XDIR]*dx;
                y0 = plo[YDIR]+iv[YDIR]*dy;
                z0 = plo[ZDIR]+iv[ZDIR]*dz;

                for(int k=0;k<2;k++)
                {
                    for(int j=0;j<2;j++)
                    {
                        for(int i=0;i<2;i++)
                        {
                            //x = x0 + (i+dist(mt))*half*dx;
                            //y = y0 + (j+dist(mt))*half*dy;
                            //z = z0 + (k+dist(mt))*half*dz;
                            x = x0 + (i+half)*half*dx;
                            y = y0 + (j+half)*half*dy;
                            z = z0 + (k+half)*half*dz;

                            if(x>=mincoords[XDIR] and x<=maxcoords[XDIR] and 
                                    y>=mincoords[YDIR] and y<=maxcoords[YDIR] and
                                    z>=mincoords[ZDIR] and z<=maxcoords[ZDIR])
                            {
                                ParticleType p = generate_particle(x,y,z,vel,
                                                 dens,eighth*dx*dy*dz,constmodel,
                                                 E,nu,bulkmod,Gama_pres,visc);
                    
                                total_mass += p.rdata(realData::mass);
                                total_vol += p.rdata(realData::volume);
                                
                                host_particles.push_back(p);
                            }
                        }
                    } 
                }
            }
        }
        
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);

    }

    // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
    // we do need this to move particles from tile 0 to the correct tile.
    Redistribute();
}

MPMParticleContainer::ParticleType MPMParticleContainer::generate_particle
        (Real x,Real y,Real z,
        Real vel[AMREX_SPACEDIM],
        Real dens, Real vol, int constmodel, Real E, Real nu,
        Real bulkmod, Real Gama_pres,Real visc)
{
    ParticleType p;
    p.id()  = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();                

    p.pos(XDIR) = x;
    p.pos(YDIR) = y;
    p.pos(ZDIR) = z;

    p.idata(intData::phase) = 0;
    p.rdata(realData::radius) = std::pow(three*fourth*vol/PI,0.33333333);

    p.rdata(realData::density) = dens;
    p.rdata(realData::xvel) = vel[XDIR];
    p.rdata(realData::yvel) = vel[YDIR];
    p.rdata(realData::zvel) = vel[ZDIR];

    p.idata(intData::constitutive_model)=constmodel;

    p.rdata(realData::E)=E;
    p.rdata(realData::nu)=nu;
    p.rdata(realData::Bulk_modulus)=bulkmod;
    p.rdata(realData::Gama_pressure)=Gama_pres;
    p.rdata(realData::Dynamic_viscosity)=visc;

    p.rdata(realData::volume)=vol;	
    p.rdata(realData::mass)=dens*vol;
    p.rdata(realData::jacobian)=1.0;
    p.rdata(realData::pressure)=0.0;
    p.rdata(realData::vol_init)=0.0;
    
    for(int comp=0;comp<NCOMP_TENSOR;comp++)
    {
        p.rdata(realData::strainrate+comp) = zero;
        p.rdata(realData::strain+comp)     = zero;
        p.rdata(realData::stress+comp)     = zero;
    }

    return(p);
}


void MPMParticleContainer::removeParticlesInsideEB()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    int lsref=mpm_ebtools::ls_refinement;

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        int np = aos.numRealParticles();

        ParticleType* pstruct = aos().dataPtr();

        amrex::Array4<amrex::Real> lsetarr=mpm_ebtools::lsphi->array(mfi);

        amrex::ParallelFor(np,[=]
                           AMREX_GPU_DEVICE (int i) noexcept
                           {
                               ParticleType& p = pstruct[i];
                               amrex::Real xp[AMREX_SPACEDIM]={p.pos(XDIR),p.pos(YDIR),p.pos(ZDIR)};

                               amrex::Real lsval=get_levelset_value(lsetarr,plo,dx,xp,lsref);

                               if(lsval<TINYVAL)
                               {
                                   p.id()=-1;
                               }

                           });
    }
    Redistribute();
}
