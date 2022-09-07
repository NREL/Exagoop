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

amrex::Real MPMParticleContainer::CalculateExactVelocity(int modenumber,amrex::Real E, amrex::Real rho, amrex::Real v0,amrex::Real L, amrex::Real time)
{
	const amrex::Real pi = 4.0*atan(1.0);
	Real beta_n = (2*modenumber-1.0)/2*pi/L;
	Real w_n = sqrt(E/rho)*beta_n;
	amrex::Real Vmex = v0/(beta_n*L)*cos(w_n*time);
	return(Vmex);
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


void MPMParticleContainer::WriteDeflectionCantilever()
{
	//Works only for serial runs

	const int lev = 0;
	auto& plev  = GetParticles(lev);

	for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
	{
		const amrex::Box& box = mfi.tilebox();
		Box nodalbox = convert(box, {1, 1, 1});

		int gid = mfi.index();
		int tid = mfi.LocalTileIndex();
		auto index = std::make_pair(gid, tid);
		auto& ptile = plev[index];
		auto& aos   = ptile.GetArrayOfStructs();
		int np = aos.numRealParticles();
		int ng =aos.numNeighborParticles();
		int nt = np+ng;

		ParticleType* pstruct = aos().dataPtr();
		amrex::ParallelFor(nt,[=]
			AMREX_GPU_DEVICE (int i) noexcept
	        {
	            ParticleType& p = pstruct[i];

	            amrex::Real xp[AMREX_SPACEDIM];

	            xp[XDIR]=p.pos(XDIR);
	            xp[YDIR]=p.pos(YDIR);
	            PrintToFile("CantileverDeflection.out")<<xp[XDIR]<<"\t"<<xp[YDIR]<<"\n";
	        });
	    }

}

void MPMParticleContainer::CalculateVelocityCantilever(Real &Vcm)
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
        return(p.rdata(realData::mass)*p.rdata(realData::yvel));
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
