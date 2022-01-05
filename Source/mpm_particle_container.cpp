#include <mpm_particle_container.H>
#include <particle_basis.H>

using namespace amrex;

void MPMParticleContainer::deposit_onto_grid(MultiFab& nodaldata) 
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();
    
    int ncomp=nodaldata.nComp();
    nodaldata.setVal(0.0);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        Box nodalbox = convert(box, {1, 1, 1});
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const int np = GetParticles(lev)[index].numRealParticles();

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
                {
                ParticleType& p = pstruct[i];
                int i_mesh,j_mesh,k_mesh;

                amrex::Real xi[AMREX_SPACEDIM];
                amrex::Real xp[AMREX_SPACEDIM];
                amrex::Real hatsize[AMREX_SPACEDIM];

                xp[XDIR]=p.pos(XDIR);
                xp[YDIR]=p.pos(YDIR);
                xp[ZDIR]=p.pos(ZDIR);

                auto iv = getParticleCell(p, plo, dxi, domain);

                hatsize[XDIR]=two*p.rdata(realData::radius);
                hatsize[YDIR]=two*p.rdata(realData::radius);
                hatsize[ZDIR]=two*p.rdata(realData::radius);
               
                for(int n=0;n<2;n++)
                {
                    for(int m=0;m<2;m++)
                    {
                        for(int l=0;l<2;l++)
                        {
                            xi[XDIR]=plo[XDIR]+(iv[0]+l)*dx[XDIR];
                            xi[YDIR]=plo[YDIR]+(iv[1]+m)*dx[YDIR];
                            xi[ZDIR]=plo[ZDIR]+(iv[2]+n)*dx[ZDIR];

                            amrex::Real basisval=hat3d(xi,xp,hatsize);
                            amrex::Real mass_contrib=p.rdata(realData::mass)*basisval*p.rdata(realData::volume);

                            amrex::Real px_contrib = p.rdata(realData::mass)*p.rdata(realData::xvel)*basisval*p.rdata(realData::volume);
                            amrex::Real py_contrib = p.rdata(realData::mass)*p.rdata(realData::yvel)*basisval*p.rdata(realData::volume);
                            amrex::Real pz_contrib = p.rdata(realData::mass)*p.rdata(realData::zvel)*basisval*p.rdata(realData::volume);

                            amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(iv[0]+l,iv[1]+m,iv[2]+n,MASS_INDEX), mass_contrib);

                            amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(iv[0]+l,iv[1]+m,iv[2]+n,VELX_INDEX), px_contrib);
                            amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(iv[0]+l,iv[1]+m,iv[2]+n,VELY_INDEX), py_contrib);
                            amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(iv[0]+l,iv[1]+m,iv[2]+n,VELZ_INDEX), pz_contrib);

                        }
                    }
                }
                });
                
               amrex::ParallelFor(
                nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                if(nodal_data_arr(i,j,k,MASS_INDEX) > 0.0)
                {
                nodal_data_arr(i,j,k,VELX_INDEX)/=nodal_data_arr(i,j,k,MASS_INDEX);
                nodal_data_arr(i,j,k,VELY_INDEX)/=nodal_data_arr(i,j,k,MASS_INDEX);
                nodal_data_arr(i,j,k,VELZ_INDEX)/=nodal_data_arr(i,j,k,MASS_INDEX);
                }


                });

    }
}

void MPMParticleContainer::moveParticles(const amrex::Real& dt,Array<Real,AMREX_SPACEDIM> gravity)
{
    BL_PROFILE("MPMParticleContainer::moveParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
    auto& plev  = GetParticles(lev);
    Real grav[]={AMREX_D_DECL(gravity[0],gravity[1],gravity[2])};

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        // now we move the particles
        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
                {
                ParticleType& p = pstruct[i];
                if(p.idata(intData::phase)!=NPOINT)
                {

                p.rdata(realData::xvel) += (p.rdata(realData::fx)/p.rdata(realData::mass) + grav[XDIR]) * dt;
                p.rdata(realData::yvel) += (p.rdata(realData::fy)/p.rdata(realData::mass) + grav[YDIR]) * dt;
                p.rdata(realData::zvel) += (p.rdata(realData::fz)/p.rdata(realData::mass) + grav[ZDIR]) * dt;

                p.pos(0) += p.rdata(realData::xvel) * dt;
                p.pos(1) += p.rdata(realData::yvel) * dt;
                p.pos(2) += p.rdata(realData::zvel) * dt;
                }
                });
    }
}

void MPMParticleContainer::writeParticles(const int n)
{
    BL_PROFILE("MPMParticleContainer::writeParticles");
    const std::string& pltfile = amrex::Concatenate("plt", n, 5);

    Vector<int> writeflags_real(realData::count,0);
    Vector<int> writeflags_int(intData::count,0);

    Vector<std::string> real_data_names;
    Vector<std::string>  int_data_names;

    real_data_names.push_back("radius");
    real_data_names.push_back("xvel");
    real_data_names.push_back("yvel");
    real_data_names.push_back("zvel");
    real_data_names.push_back("fx");
    real_data_names.push_back("fy");
    real_data_names.push_back("fz");
    real_data_names.push_back("volume");
    real_data_names.push_back("mass");
    real_data_names.push_back("density");
    int_data_names.push_back("phase");

    writeflags_int[intData::phase]=1;

    writeflags_real[realData::radius]=1;
    writeflags_real[realData::xvel]=1;
    writeflags_real[realData::yvel]=1;
    writeflags_real[realData::zvel]=1;
    writeflags_real[realData::fx]=1;
    writeflags_real[realData::fy]=1;
    writeflags_real[realData::fz]=1;
    writeflags_real[realData::mass]=1;

    WritePlotFile(pltfile, "particles",writeflags_real, 
            writeflags_int, real_data_names, int_data_names);
}
