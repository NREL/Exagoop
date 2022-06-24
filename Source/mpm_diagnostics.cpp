#include <mpm_particle_container.H>
#include <interpolants.H>

void MPMParticleContainer::CalculateEnergies(Real &TKE,Real &TSE)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    TKE=0.0;
    TSE=0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    TKE = amrex::ReduceSum(*this, [=] 
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real 
    {
        return(0.5*p.rdata(realData::mass)*
               (p.rdata(realData::xvel)*p.rdata(realData::xvel)+
                p.rdata(realData::yvel)*p.rdata(realData::yvel)+
                p.rdata(realData::zvel)*p.rdata(realData::zvel)) );
    });

    TSE = amrex::ReduceSum(*this, [=] 
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real 
    {
        return(0.5*p.rdata(realData::volume)*
               (p.rdata(realData::stress+XX)*p.rdata(realData::strain+XX)+
                p.rdata(realData::stress+YY)*p.rdata(realData::strain+YY)+
                p.rdata(realData::stress+ZZ)*p.rdata(realData::strain+ZZ)+
                p.rdata(realData::stress+XY)*p.rdata(realData::strain+XY)*2.0+
                p.rdata(realData::stress+YZ)*p.rdata(realData::strain+YZ)*2.0+
                p.rdata(realData::stress+XZ)*p.rdata(realData::strain+XZ)*2.0));
    });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealSum(TKE);
    ParallelDescriptor::ReduceRealSum(TSE);
#endif

}

void MPMParticleContainer::CalculateVelocity(Real &Vcm)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    Real Vcmx=0.0;
    Real mass_tot=0.0;
    
    using PType = typename MPMParticleContainer::SuperParticleType;
    Vcmx = amrex::ReduceSum(*this, [=] 
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real 
    {
        return(p.rdata(realData::mass)*p.rdata(realData::xvel));
    });
    
    mass_tot = amrex::ReduceSum(*this, [=] 
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real 
    {
        return(p.rdata(realData::mass));
    });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealSum(Vcmx);
    ParallelDescriptor::ReduceRealSum(mass_tot);
#endif
    
    Vcm=Vcmx/mass_tot;
}

void MPMParticleContainer::FindWaterFront(Real &Xwf)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    Real wf_x=0.0;
    Real mass_tot=0.0;
    
    using PType = typename MPMParticleContainer::SuperParticleType;
    wf_x = amrex::ReduceMax(*this, [=] 
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real 
    {
        return(p.pos(XDIR));
    });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealMax(wf_x);
#endif

    Xwf=wf_x;
}
