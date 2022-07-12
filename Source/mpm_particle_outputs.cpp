#include <mpm_particle_container.H>
#include <interpolants.H>

void MPMParticleContainer::update_density_field(MultiFab& nodaldata,int refratio,Real smoothfactor)
{
    nodaldata.setVal(zero);
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
        Box nodalbox = convert(refbox, {1, 1, 1});
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        int np = aos.numRealParticles();
        int ng = aos.numNeighborParticles();
        int nt = np+ng;

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(nt,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            auto iv = getParticleCell(p, plo, dxi, domain);

            for(int n=0;n<2;n++)
            {
                for(int m=0;m<2;m++)
                {
                    for(int l=0;l<2;l++)
                    {
                        IntVect ivlocal(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n);

                        if(nodalbox.contains(ivlocal))
                        {
                            amrex::Real xp[AMREX_SPACEDIM];
                            amrex::Real xi[AMREX_SPACEDIM];
                            amrex::Real weight;

                            xp[XDIR]=p.pos(XDIR);
                            xp[YDIR]=p.pos(YDIR);
                            xp[ZDIR]=p.pos(ZDIR);

                            xi[XDIR]=plo[XDIR]+ivlocal[XDIR]*dx[XDIR];
                            xi[YDIR]=plo[YDIR]+ivlocal[YDIR]*dx[YDIR];
                            xi[ZDIR]=plo[ZDIR]+ivlocal[ZDIR]*dx[ZDIR];

                            weight=p.rdata(realData::mass)*spherical_gaussian(xi,xp,smoothfactor*p.rdata(realData::radius));

                            amrex::Gpu::Atomic::AddNoRet(
                                &nodal_data_arr(ivlocal),
                                weight);
                        }
                    }
                }
            }

        });

    }
}

void MPMParticleContainer::writeParticles(const int n)
{
    BL_PROFILE("MPMParticleContainer::writeParticles");
    const std::string& pltfile = amrex::Concatenate("plt", n, 5);

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
    real_data_names.push_back("Bulk_modulous");
    real_data_names.push_back("Gama_pressure");
    real_data_names.push_back("Dynamic_viscosity");
    real_data_names.push_back("void_ratio");

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
    writeflags_real[realData::Bulk_modulous]=0;
    writeflags_real[realData::Gama_pressure]=0;
    writeflags_real[realData::Dynamic_viscosity]=0;
    writeflags_real[realData::void_ratio]=1;
    
    WritePlotFile(pltfile, "particles",writeflags_real, 
                  writeflags_int, real_data_names, int_data_names);
}
