#include <mpm_particle_container.H>
#include <particle_basis.H>

using namespace amrex;
using namespace std;

void MPMParticleContainer::deposit_onto_nodes() 
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    //zero data
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
        {
                ParticleType& p = pstruct[i];

                if(p.idata(intData::phase)==NPOINT)
                {
                   p.rdata(realData::fx) = zero;
                   p.rdata(realData::fy) = zero;
                   p.rdata(realData::fz) = zero;

                   p.rdata(realData::xvel) = zero;
                   p.rdata(realData::xvel) = zero;
                   p.rdata(realData::xvel) = zero;

                   p.rdata(realData::mass)    = zero;
                   p.rdata(realData::volume)  = zero;
                   p.rdata(realData::density) = zero;
                }
        });
    }

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p1 = pstruct[i];

            if(p1.idata(intData::phase)==NPOINT)
            {
                for (const auto& p2 : nbor_data.getNeighbors(i))
                { 
                    if(p2.idata(intData::phase)!=NPOINT)
                    {
                        amrex::Real xi[AMREX_SPACEDIM];
                        amrex::Real xp[AMREX_SPACEDIM];
                        amrex::Real hatsize[AMREX_SPACEDIM];

                        xi[XDIR]=p1.pos(XDIR);
                        xi[YDIR]=p1.pos(YDIR);
                        xi[ZDIR]=p1.pos(ZDIR);

                        xp[XDIR]=p2.pos(XDIR);
                        xp[YDIR]=p2.pos(YDIR);
                        xp[ZDIR]=p2.pos(ZDIR);

                        hatsize[XDIR]=two*p2.rdata(realData::radius);
                        hatsize[YDIR]=two*p2.rdata(realData::radius);
                        hatsize[ZDIR]=two*p2.rdata(realData::radius);

                        amrex::Real basisval=hat3d(xi,xp,hatsize);
                        p1.rdata(realData::mass) += p2.rdata(realData::mass)*basisval*p2.rdata(realData::volume);
                        p1.rdata(realData::xvel) += p2.rdata(realData::mass)*p2.rdata(realData::xvel)*basisval*p2.rdata(realData::volume);
                        p1.rdata(realData::yvel) += p2.rdata(realData::mass)*p2.rdata(realData::yvel)*basisval*p2.rdata(realData::volume);
                        p1.rdata(realData::zvel) += p2.rdata(realData::mass)*p2.rdata(realData::zvel)*basisval*p2.rdata(realData::volume);
                    }
                }        
            }
        });
        
        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];

            if(p.idata(intData::phase)==NPOINT)
            {
                if(p.rdata(realData::mass)!=0)
                {
                   p.rdata(realData::xvel)=p.rdata(realData::xvel)/p.rdata(realData::mass);
                   p.rdata(realData::yvel)=p.rdata(realData::yvel)/p.rdata(realData::mass);
                   p.rdata(realData::zvel)=p.rdata(realData::zvel)/p.rdata(realData::mass);
                }
            }
        });
    }
}

void MPMParticleContainer::deposit_onto_grid(MultiFab& nodaldata) 
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    //zero data
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        Box nodalbox = convert(box, {1, 1, 1});
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
                {
                ParticleType& p = pstruct[i];
                int i_mesh,j_mesh,k_mesh;
                amrex::Real tol=TINYVAL;

                amrex::Real xi[AMREX_SPACEDIM];
                amrex::Real xp[AMREX_SPACEDIM];
                amrex::Real hatsize[AMREX_SPACEDIM];

                xp[XDIR]=p.pos(XDIR);
                xp[YDIR]=p.pos(YDIR);
                xp[ZDIR]=p.pos(ZDIR);

                i_mesh=amrex::Math::floor((p.pos(XDIR)-plo[XDIR]+tol)/dx[XDIR]);
                j_mesh=amrex::Math::floor((p.pos(YDIR)-plo[YDIR]+tol)/dx[YDIR]);
                k_mesh=amrex::Math::floor((p.pos(ZDIR)-plo[ZDIR]+tol)/dx[ZDIR]);

                hatsize[XDIR]=two*p.rdata(realData::radius);
                hatsize[YDIR]=two*p.rdata(realData::radius);
                hatsize[ZDIR]=two*p.rdata(realData::radius);

                for(int n=0;n<2;n++)
                {
                    for(int m=0;m<2;m++)
                    {
                        for(int l=0;l<2;l++)
                        {
                            xi[XDIR]=plo[XDIR]+(i_mesh+l)*dx[XDIR];
                            xi[YDIR]=plo[YDIR]+(j_mesh+m)*dx[YDIR];
                            xi[ZDIR]=plo[ZDIR]+(k_mesh+n)*dx[ZDIR];

                            amrex::Real basisval=hat3d(xi,xp,hatsize);
                            basisval=1.0;
                            amrex::Real mass_contrib=p.rdata(realData::mass)*basisval*p.rdata(realData::volume);

                            amrex::Real px_contrib = p.rdata(realData::mass)*p.rdata(realData::xvel)*basisval*p.rdata(realData::volume);
                            amrex::Real py_contrib = p.rdata(realData::mass)*p.rdata(realData::yvel)*basisval*p.rdata(realData::volume);
                            amrex::Real pz_contrib = p.rdata(realData::mass)*p.rdata(realData::zvel)*basisval*p.rdata(realData::volume);

                            amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(i_mesh,j_mesh,k_mesh,MASS_INDEX), mass_contrib);

                            amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(i_mesh,j_mesh,k_mesh,VELX_INDEX), px_contrib);
                            amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(i_mesh,j_mesh,k_mesh,VELY_INDEX), py_contrib);
                            amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(i_mesh,j_mesh,k_mesh,VELZ_INDEX), pz_contrib);

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
