#include <mpm_particle_container.H>

using namespace amrex;
using namespace std;

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

            p.rdata(realData::xvel) += (p.rdata(realData::fx)/p.rdata(realData::mass) + grav[XDIR]) * dt;
            p.rdata(realData::yvel) += (p.rdata(realData::fy)/p.rdata(realData::mass) + grav[YDIR]) * dt;
            p.rdata(realData::zvel) += (p.rdata(realData::fz)/p.rdata(realData::mass) + grav[ZDIR]) * dt;

            p.pos(0) += p.rdata(realData::xvel) * dt;
            p.pos(1) += p.rdata(realData::yvel) * dt;
            p.pos(2) += p.rdata(realData::zvel) * dt;
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

    writeflags_real[realData::radius]=1;
    writeflags_real[realData::xvel]=1;
    writeflags_real[realData::yvel]=1;
    writeflags_real[realData::zvel]=1;
    writeflags_real[realData::fx]=1;
    writeflags_real[realData::fy]=1;
    writeflags_real[realData::fz]=1;

    WritePlotFile(pltfile, "particles",writeflags_real, 
            writeflags_int, real_data_names, int_data_names);
}
