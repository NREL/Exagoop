#include <mpm_particle_container.H>
#include <interpolants.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_AmrMesh.H>

void MPMParticleContainer::update_density_field(MultiFab& densdata,int refratio,Real smoothfactor)
{
    int ng_dens=3;
    densdata.setVal(zero,ng_dens);
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    GpuArray<Real,AMREX_SPACEDIM> dxi = geom.InvCellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    Box domain = geom.Domain();
    domain.refine(refratio);

    dxi[XDIR]*=refratio;
    dxi[YDIR]*=refratio;
    dxi[ZDIR]*=refratio;

    dx[XDIR]/=refratio;
    dx[YDIR]/=refratio;
    dx[ZDIR]/=refratio;


    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        amrex::Box box = mfi.tilebox();
        amrex::Box& refbox = box.refine(refratio);
        const amrex::Box& refboxgrow = amrex::grow(refbox,ng_dens);
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        int np = aos.numRealParticles();
        int nt = np;

        Array4<Real> dens_data_arr=densdata.array(mfi);

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(nt,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            auto iv = getParticleCell(p, plo, dxi, domain);

            for(int n=-3;n<=3;n++)
            {
                for(int m=-3;m<=3;m++)
                {
                    for(int l=-3;l<=3;l++)
                    {
                        IntVect ivlocal(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n);

                        if(refboxgrow.contains(ivlocal))
                        {
                            amrex::Real xp[AMREX_SPACEDIM];
                            amrex::Real xi[AMREX_SPACEDIM];
                            amrex::Real weight;

                            xp[XDIR]=p.pos(XDIR);
                            xp[YDIR]=p.pos(YDIR);
                            xp[ZDIR]=p.pos(ZDIR);

                            xi[XDIR]=plo[XDIR]+(ivlocal[XDIR]+half)*dx[XDIR];
                            xi[YDIR]=plo[YDIR]+(ivlocal[YDIR]+half)*dx[YDIR];
                            xi[ZDIR]=plo[ZDIR]+(ivlocal[ZDIR]+half)*dx[ZDIR];

                            weight=p.rdata(realData::mass)*
                            spherical_gaussian(xi,xp,smoothfactor*p.rdata(realData::radius));

                            amrex::Gpu::Atomic::AddNoRet(
                                &dens_data_arr(ivlocal),
                                weight);
                        }
                        //else
                        //{
                        //   amrex::Print()<<"iv,box,p:"<<iv<<"\t"<<box<<"\t"<<p.pos(0)<<"\t"<<p.pos(1)<<"\t"<<p.pos(2)<<"\n";
                        //}
                    }
                }
            }

        });

    }

    densdata.SumBoundary(geom.periodicity());
}

void MPMParticleContainer::writeParticles(std::string prefix_particlefilename, int num_of_digits_in_filenames, const int n)
{
    BL_PROFILE("MPMParticleContainer::writeParticles");
    const std::string& pltfile = amrex::Concatenate(prefix_particlefilename, n, num_of_digits_in_filenames);

    Vector<int> writeflags_real(realData::count,1);
    Vector<int> writeflags_int(intData::count,0);


    Vector<std::string> real_data_names;
    Vector<std::string>  int_data_names;

    real_data_names.push_back("radius");
    real_data_names.push_back("xvel");
    real_data_names.push_back("yvel");
    real_data_names.push_back("zvel");
    for(int i=0;i<6;i++)
    {
        real_data_names.push_back(amrex::Concatenate("strainrate_", i, 1));
    }
    for(int i=0;i<6;i++)
    {
        real_data_names.push_back(amrex::Concatenate("strain_", i, 1));
    }
    for(int i=0;i<6;i++)
    {
        real_data_names.push_back(amrex::Concatenate("stress_", i, 1));
    }
    real_data_names.push_back("volume");
    real_data_names.push_back("mass");
    real_data_names.push_back("density");
    real_data_names.push_back("jacobian");
    real_data_names.push_back("pressure");
    real_data_names.push_back("vol_init");
    real_data_names.push_back("E");
    real_data_names.push_back("nu");
    real_data_names.push_back("Bulk_modulus");
    real_data_names.push_back("Gama_pressure");
    real_data_names.push_back("Dynamic_viscosity");

    int_data_names.push_back("phase");
    int_data_names.push_back("constitutive_model");


    writeflags_int[intData::phase]=1;
    writeflags_int[intData::constitutive_model]=1;

    writeflags_real[realData::radius]=1;
    writeflags_real[realData::xvel]=1;
    writeflags_real[realData::yvel]=1;
    writeflags_real[realData::zvel]=1;
    writeflags_real[realData::mass]=1;
    writeflags_real[realData::jacobian]=1;
    writeflags_real[realData::pressure]=1;
    writeflags_real[realData::vol_init]=1;
    writeflags_real[realData::E]=0;
    writeflags_real[realData::nu]=0;
    writeflags_real[realData::Bulk_modulus]=0;
    writeflags_real[realData::Gama_pressure]=0;
    writeflags_real[realData::Dynamic_viscosity]=0;
    
    WritePlotFile(pltfile, "particles",writeflags_real, 
                  writeflags_int, real_data_names, int_data_names);
}

void MPMParticleContainer::WriteHeader(const std::string& name, bool is_checkpoint, amrex::Real cur_time, int nstep, int EB_generate_max_level, int output_it) const
{
    if(ParallelDescriptor::IOProcessor())
    {
    	const int finest_level = 0;
        std::string HeaderFileName(name + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;

        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        HeaderFile.open(HeaderFileName.c_str(),
                        std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

        if(!HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);
        if(is_checkpoint) {
            HeaderFile << "Checkpoint version: 1\n";
        } else {
            HeaderFile << "HyperCLaw-V1.1\n";
        }

        HeaderFile << nstep << "\n";
        HeaderFile << output_it << "\n";

#ifdef AMREX_USE_EB
        HeaderFile << EB_generate_max_level << "\n";
#endif

        HeaderFile << cur_time << "\n";


    }
}

void MPMParticleContainer::writeCheckpointFile(std::string prefix_particlefilename, int num_of_digits_in_filenames, amrex::Real cur_time, int nstep, int output_it)
{
	BL_PROFILE("MPMParticleContainer::writeCheckpointFile");
	const int m_nstep = output_it;
	const int finest_level=0;
	const int EB_generate_max_level=0;
	std::string level_prefix = "Level_";
	const std::string& checkpointname = amrex::Concatenate(prefix_particlefilename, m_nstep, num_of_digits_in_filenames);
	amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level + 1, true);
	bool is_checkpoint = true;
	WriteHeader(checkpointname, is_checkpoint, cur_time, nstep, EB_generate_max_level,output_it);

	amrex::Vector<std::string> real_data_names;
	real_data_names.push_back("radius");
	real_data_names.push_back("xvel");
	real_data_names.push_back("yvel");
	real_data_names.push_back("zvel");
	for(int i=0;i<6;i++)
	{
		real_data_names.push_back(amrex::Concatenate("strainrate_", i, 1));
	}
	for(int i=0;i<6;i++)
	{
		real_data_names.push_back(amrex::Concatenate("strain_", i, 1));
	}
	for(int i=0;i<6;i++)
	{
		real_data_names.push_back(amrex::Concatenate("stress_", i, 1));
	}
	real_data_names.push_back("volume");
	real_data_names.push_back("mass");
	real_data_names.push_back("density");
	real_data_names.push_back("jacobian");
	real_data_names.push_back("pressure");
	real_data_names.push_back("vol_init");
	real_data_names.push_back("E");
	real_data_names.push_back("nu");
	real_data_names.push_back("Bulk_modulus");
	real_data_names.push_back("Gama_pressure");
	real_data_names.push_back("Dynamic_viscosity");

	amrex::Vector<std::string> int_data_names;
	int_data_names.push_back("phase");
	int_data_names.push_back("constitutive_model");

	Checkpoint( checkpointname, "particles", is_checkpoint, real_data_names, int_data_names);
}

void GotoNextLine(std::istream& is)
{
       constexpr std::streamsize bl_ignore_max{100000};
           is.ignore(bl_ignore_max, '\n');
}

void MPMParticleContainer::readCheckpointFile(std::string & restart_chkfile, int &nstep, double &cur_time, int &output_it)
{
	BL_PROFILE("MPMParticleContainer::readCheckpointFile");

	amrex::Print() << "Restarting from checkpoint " << restart_chkfile << "\n";

	Real prob_lo[AMREX_SPACEDIM];
	Real prob_hi[AMREX_SPACEDIM];
	const int max_level = 0;
	const int finest_level=0;

	/***************************************************************************
	   ** Load header: set up problem domain (including BoxArray)                 *
	   **              allocate PeleLM memory (PeleLM::AllocateArrays)            *
	   **              (by calling MakeNewLevelFromScratch)                       *
	   ****************************************************************************/

	std::string File(restart_chkfile + "/Header");

	VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

	Vector<char> fileCharPtr;
	ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
	std::string fileCharPtrString(fileCharPtr.dataPtr());
	std::istringstream is(fileCharPtrString, std::istringstream::in);

	std::string line, word;

	// Start reading from checkpoint file

	// Title line
	std::getline(is, line);

	// Finest level
	int chk_finest_level = 0;

	// Step count
	is >> nstep;
	GotoNextLine(is);

	// Output number for plot files
	is >> output_it;
	GotoNextLine(is);

	#ifdef AMREX_USE_EB
	   // Finest level at which EB was generated
	   // actually used independently, so just skip ...
	   std::getline(is, line);

	   // ... but to be backward compatible, if we get a float,
	   // let's assume it's m_cur_time
	   if (line.find('.') != std::string::npos) {
	      cur_time = std::stod(line);

	   } else {
	      // Skip line and read current time
	      is >> cur_time;
	      GotoNextLine(is);
	   }
	#else

	   // Current time
	   is >> cur_time;
	   GotoNextLine(is);
	#endif

	   Restart(restart_chkfile,"particles", true);

	   if (m_verbose) {
	      amrex::Print() << "Restart complete" << std::endl;
	   }

}
