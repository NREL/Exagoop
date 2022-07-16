#include <mpm_particle_container.H>
#include <constants.H>

void MPMParticleContainer::InitParticles (const std::string& filename,Real *total_mass,Real *total_vol)
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

        *total_mass=0.0;
        *total_vol=0.0;

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
            ifs >> p.idata(intData::constitutive_model);		//Commented only for getting the HPRO inputfile to work
            //ifs>>junk;
            //p.idata(intData::constitutive_model)=0;
            if(p.idata(intData::constitutive_model)==0)	//Elastic solid
            {
            	ifs >> p.rdata(realData::E);
            	ifs >> p.rdata(realData::nu);
            	p.rdata(realData::Bulk_modulus)=0.0;
            	p.rdata(realData::Gama_pressure)=0.0;
            	p.rdata(realData::Dynamic_viscosity)=0.0;
                p.rdata(realData::void_ratio)=0.0;
            }
            else if(p.idata(intData::constitutive_model)==1)
            {
            	p.rdata(realData::E)=0.0;
            	p.rdata(realData::nu)=0.0;
            	ifs >> p.rdata(realData::Bulk_modulus);
            	ifs >> p.rdata(realData::Gama_pressure);
            	ifs >> p.rdata(realData::Dynamic_viscosity);
                p.rdata(realData::void_ratio)=0.0;
            }
            else if(p.idata(intData::constitutive_model)==2) // Yudong: hypoplastic model
            {
            	p.rdata(realData::E)=0.0;
            	p.rdata(realData::nu)=0.0;
                p.rdata(realData::Bulk_modulus)=0.0;
            	p.rdata(realData::Gama_pressure)=0.0;
            	p.rdata(realData::Dynamic_viscosity)=0.0;
            	ifs >> p.rdata(realData::void_ratio);
            }
            else
            {
            	amrex::Abort("\nIncorrect constitutive model. Please check your particle file");
            }

            // Set other particle properties
            p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three);	//This is a dummy initialisation. We will correct this value later
            p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);

            *total_mass +=p.rdata(realData::mass);
            *total_vol +=p.rdata(realData::volume);
            p.rdata(realData::jacobian)	   = 1.0;
            p.rdata(realData::vol_init)	   = 0.0;

            double initial_strainrate = zero;
            double initial_spinrate = zero;
            double initial_strain = zero;
            double initial_stress = zero;
            
            if(p.idata(intData::constitutive_model)==2){
                initial_strainrate = zero;
                initial_spinrate = zero;
                initial_strain = zero;
                initial_stress = -10; // initial stress provide by input
            }
            
            for(int comp=0;comp<NCOMP_TENSOR;comp++)
            {
                p.rdata(realData::strainrate+comp) = initial_strainrate;
                p.rdata(realData::spinrate+comp) = initial_spinrate;
                p.rdata(realData::strain+comp)     = initial_strain;
                p.rdata(realData::stress+comp)     = initial_stress; 
            }
            
            host_particles.push_back(p);

            //amrex::Print()<<"\n constitutive model = "<<p.idata(intData::constitutive_model)<<", init_stress_xx = "<<p.rdata(realData::stress);

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
