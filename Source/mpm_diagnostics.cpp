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

void MPMParticleContainer::CalculateSurfaceIntegralTop(Array<Real,AMREX_SPACEDIM> gravity, Real &Fy_top, Real &Fy_bottom)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    Real Mvy=0.0;
    Real Fg=0.0;
    Fy_bottom=0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    Mvy = amrex::ReduceSum(*this, [=]
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
    {
        return(p.rdata(realData::mass)*p.rdata(realData::yacceleration));
    });

    Fg = amrex::ReduceSum(*this, [=]
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
    {
        return(p.rdata(realData::mass)*gravity[YDIR]);
    });


#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealSum(Mvy);
    ParallelDescriptor::ReduceRealSum(Fg);
#endif

    Fy_top = Mvy+fabs(Fg)+fabs(Fy_bottom);
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

void MPMParticleContainer::CalculateErrorTVB(Real tvb_E,Real tvb_v0,Real tvb_L,Real tvb_rho,Real err)
{
	const Real pi = atan(1.0)*4.0;
	Real c = sqrt(tvb_E/tvb_rho);
	Real w0 = c*pi/tvb_L;

}

void MPMParticleContainer::CalculateErrorP2G(MultiFab& nodaldata,amrex::Real p2g_L,amrex::Real p2g_f, int ncell)
{
	const int lev = 0;
	    const Geometry& geom = Geom(lev);
	    auto& plev  = GetParticles(lev);
	    const auto dxi = geom.InvCellSizeArray();
	    const auto dx = geom.CellSizeArray();
	    const auto plo = geom.ProbLoArray();
	    const auto domain = geom.Domain();
	    std::string outputfile;

	    const int* lo = domain.loVect ();
	    const int* hi = domain.hiVect ();
	    outputfile = amrex::Concatenate("P2GTest1",ncell,3);

	    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
	    {
	        //already nodal as mfi is from nodaldata
	        const Box& nodalbox=mfi.validbox();

	        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

	        amrex::ParallelFor(nodalbox,[=]
	        AMREX_GPU_DEVICE (int i,int j,int k) noexcept
	        {
	        	amrex::Real x,v;
	        	x = lo[0]+i*dx[0];
	        	v= nodal_data_arr(i,0,0,VELY_INDEX);
	        	if(j==0 and k==0)
	        		{
	        	PrintToFile(outputfile)<<x<<"\t"<<v<<"\n";
	        		}
	        });
	    }
	    outputfile = amrex::Concatenate("P2GTest2",ncell,3);
	    //std::ofstream ofs(outputfile, std::ofstream::out);

	      /*amrex::Print(ofs)
	        << "L, rho, umax, p, T, gamma, mu, k, Re, Ma, Pr, dpdx, G, radius"
	        << std::endl;
	      amrex::Print(ofs).SetPrecision(17)
	        << L << "," << PeleC::h_prob_parm_device->rho << ","
	        << PeleC::h_prob_parm_device->umax << "," << PeleC::h_prob_parm_device->p
	        << "," << PeleC::h_prob_parm_device->T << "," << eos.gamma << ","
	        << trans_parm.const_viscosity << "," << trans_parm.const_conductivity << ","
	        << PeleC::h_prob_parm_device->Re << "," << PeleC::h_prob_parm_device->Ma
	        << "," << PeleC::h_prob_parm_device->Pr << ","
	        << PeleC::h_prob_parm_device->dpdx << "," << PeleC::h_prob_parm_device->G
	        << "," << PeleC::h_prob_parm_device->radius << std::endl;*/

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
	    		Real y_exact;

	    		amrex::Real xp[AMREX_SPACEDIM];

	    		xp[XDIR]=p.pos(XDIR);

	    		y_exact = p.rdata(realData::yvel);
	    		PrintToFile(outputfile).SetPrecision(17)<<xp[XDIR]<<"\t"<<y_exact<<"\n";
	    		//amrex::Print(ofs).SetPrecision(17)<<xp[XDIR]<<"\t"<<y_exact<<"\n";

	    	});
	    }
}

void MPMParticleContainer::WriteDeflectionTVB(Real tvb_E,Real tvb_v0,Real tvb_L,Real tvb_rho, Real time, int output_it)
{
	//Works only for serial runs

	const int lev = 0;
	auto& plev  = GetParticles(lev);
	Real c = sqrt(tvb_E/(tvb_rho));
	const Real pi = atan(1.0)*4.0;
	Real w0 = c*pi/tvb_L;
	Real Amplitude = tvb_v0*tvb_L/(c*pi)*sin(w0*time);
	std::string outputfile;
	int num_of_digits_in_filenames = 5;
	outputfile = amrex::Concatenate("TVB_Deflection_", output_it, num_of_digits_in_filenames);
	amrex::Print()<<"\nE = "<<tvb_E<<" "<<tvb_rho<<" "<<c<<" "<<w0;

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
	            Real y_exact;

	            amrex::Real xp[AMREX_SPACEDIM];

	            xp[XDIR]=p.pos(XDIR);
	            xp[YDIR]=p.pos(YDIR);
	            y_exact = Amplitude*sin(pi*xp[XDIR]/tvb_L);
	            PrintToFile(outputfile)<<xp[XDIR]<<"\t"<<xp[YDIR]<<"\t"<<y_exact<<"\n";
	        });
	    }

}

amrex::Real MPMParticleContainer::CalculateEffectiveSpringConstant(amrex::Real Area, amrex::Real L0)
{
	//First calculate the total strain energy
	const int lev = 0;
	const Geometry& geom = Geom(lev);
	auto& plev  = GetParticles(lev);
	const auto dxi = geom.InvCellSizeArray();
	const auto dx = geom.CellSizeArray();
	const auto plo = geom.ProbLoArray();
	const auto domain = geom.Domain();

	amrex::Real TSE=0.0;
	amrex::Real Total_vol=0.0;
	amrex::Real deflection = 0.0;
	amrex::Real Restoring_force = 0.0;
	amrex::Real smallval = 1e-10;
	amrex::Real Calculated_Spring_Const = 0.0;

	using PType = typename MPMParticleContainer::SuperParticleType;

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
	    ParallelDescriptor::ReduceRealSum(TSE);
	#endif

	//Then Calculate the total volume at this instant
	Total_vol = amrex::ReduceSum(*this, [=]
	    	    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
	{
		return(p.rdata(realData::volume));
	});

	//Calculate the deflection
	deflection = L0-Total_vol/Area;

	//Calculate the spring constant
	if(fabs(deflection)<=smallval)
	{
		Restoring_force=0.0;
	}
	else
	{
		Restoring_force= 2*TSE/deflection;
		Calculated_Spring_Const = 2*TSE/(deflection*deflection);

	}

	PrintToFile("SpringConst.out")<<Calculated_Spring_Const<<"\n";

	//Calculate and return the restoring force
	return(Restoring_force);

}


