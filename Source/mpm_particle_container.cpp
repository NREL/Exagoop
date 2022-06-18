#include <mpm_particle_container.H>
#include <interpolants.H>
#include <constitutive_models.H>

using namespace amrex;

void MPMParticleContainer::apply_constitutive_model(const amrex::Real& dt,
                                                    amrex::Real applied_strainrate=0.0
													)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        int np = aos.numRealParticles();
        int ng = aos.numNeighborParticles();
        int nt = np+ng;

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(nt,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];

            amrex::Real xp[AMREX_SPACEDIM];
            amrex::Real strainrate[NCOMP_TENSOR];
            amrex::Real strain[NCOMP_TENSOR];
            amrex::Real stress[NCOMP_TENSOR];

            for(int d=0;d<NCOMP_TENSOR;d++)
            {
                p.rdata(realData::strain+d) += dt*p.rdata(realData::strainrate+d);
            }
            //apply axial strain
            p.rdata(realData::strain+XX) += dt*applied_strainrate;
            p.rdata(realData::strain+YY) += dt*applied_strainrate;
            p.rdata(realData::strain+ZZ) += dt*applied_strainrate;

            for(int d=0;d<NCOMP_TENSOR;d++)
            {
                strainrate[d]=p.rdata(realData::strainrate+d);
                strain[d]=p.rdata(realData::strain+d);
            }


            if(p.idata(intData::constitutive_model)==0)		//Elastic solid
            {
            	linear_elastic(strain,strainrate,stress,p.rdata(realData::E),p.rdata(realData::nu));
            }
            else if(p.idata(intData::constitutive_model==1))		//Viscous fluid with approximate EoS
            {
            	p.rdata(realData::pressure) = p.rdata(realData::Bulk_modulous)*(pow(1/p.rdata(realData::jacobian),p.rdata(realData::Gama_pressure))-1.0);
            	Newtonian_Fluid(strainrate,stress,p.rdata(realData::Dynamic_viscosity),p.rdata(realData::pressure));
            }

            for(int d=0;d<NCOMP_TENSOR;d++)
            {
                p.rdata(realData::stress+d)=stress[d];
            }
        });
    }
} 

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
		amrex::ParallelFor(nt,[=,&TKE,&TSE]
			AMREX_GPU_DEVICE (int i) noexcept
	        {
	            ParticleType& p = pstruct[i];
	            TKE += 0.5*p.rdata(realData::mass)*(p.rdata(realData::xvel)*p.rdata(realData::xvel)+
	            		p.rdata(realData::yvel)*p.rdata(realData::yvel)+
						p.rdata(realData::zvel)*p.rdata(realData::zvel));
	            TSE += 0.5*p.rdata(realData::volume)*
	            		(p.rdata(realData::stress+XX)*p.rdata(realData::strain+XX)+
	            		 p.rdata(realData::stress+YY)*p.rdata(realData::strain+YY)+
						 p.rdata(realData::stress+ZZ)*p.rdata(realData::strain+ZZ)+
						 p.rdata(realData::stress+XY)*p.rdata(realData::strain+XY)*2.0+
						 p.rdata(realData::stress+YZ)*p.rdata(realData::strain+YZ)*2.0+
						 p.rdata(realData::stress+XZ)*p.rdata(realData::strain+XZ)*2.0);

	        });
	}

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

	Vcm=0.0;
	Real mass_tot=0.0;



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
		amrex::ParallelFor(nt,[=,&Vcm,&mass_tot]
			AMREX_GPU_DEVICE (int i) noexcept
	        {
	            ParticleType& p = pstruct[i];
	            Vcm += p.rdata(realData::mass)*p.rdata(realData::yvel);
	            mass_tot +=p.rdata(realData::mass);

	        });
	}
	Vcm=Vcm/mass_tot;

}

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


amrex::Real MPMParticleContainer::Calculate_time_step()
{
	const int lev = 0;
	const Geometry& geom = Geom(lev);
	auto& plev  = GetParticles(lev);
	const auto dx = geom.CellSizeArray();
	amrex::Real dt = std::numeric_limits<amrex::Real>::max();


	for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
	{
		int gid = mfi.index();
		int tid = mfi.LocalTileIndex();
		auto index = std::make_pair(gid, tid);

		auto& ptile = plev[index];
		auto& aos   = ptile.GetArrayOfStructs();
		int np = aos.numRealParticles();
		int ng =aos.numNeighborParticles();
		int nt = np+ng;

		ParticleType* pstruct = aos().dataPtr();
		amrex::ParallelFor(nt,[=, &dt]
	        AMREX_GPU_DEVICE (int i) noexcept
	        {
			ParticleType& p = pstruct[i];
			amrex::Real Cs;
			amrex::Real lambda;
			amrex::Real mu;


			if(p.idata(intData::constitutive_model)==1)
			{
				Cs = sqrt(p.rdata(realData::Bulk_modulous)/p.rdata(realData::density));
			}
			else if(p.idata(intData::constitutive_model)==0)
			{
				lambda=p.rdata(realData::E)*p.rdata(realData::nu)/((1+p.rdata(realData::nu))*(1-2.0*p.rdata(realData::nu)));
				mu=p.rdata(realData::E)/(2.0*(1+p.rdata(realData::nu)));
				Cs = sqrt((lambda+2.0*mu)/p.rdata(realData::density));
			}

			//amrex::Print()<<"\n "<<p.rdata(realData::xvel)<<" "<<p.rdata(realData::yvel)<<" "<<p.rdata(realData::zvel)<<" "<<Cs;

			AMREX_D_TERM(const amrex::Real dt1 = dx[0] / (Cs + amrex::Math::abs(p.rdata(realData::xvel)));
			                   dt = amrex::min<amrex::Real>(dt, dt1);,
			             const amrex::Real dt2 = dx[1] / (Cs + amrex::Math::abs(p.rdata(realData::yvel)));
			                   dt = amrex::min<amrex::Real>(dt, dt2);,
			             const amrex::Real dt3 = dx[2] / (Cs + amrex::Math::abs(p.rdata(realData::zvel)));
			                   dt = amrex::min<amrex::Real>(dt, dt3););

	        });
	}
	if(dt<1e-10)
	{
		amrex::Print()<<"\nWarning: Time step is getting too low (dt = "<<dt<<" )";
	}
	return(dt);
}

void MPMParticleContainer::deposit_onto_grid(MultiFab& nodaldata,
                                             Array<Real,AMREX_SPACEDIM> gravity,
                                             int external_loads_present,
                                             Array<Real,AMREX_SPACEDIM> force_slab_lo,
                                             Array<Real,AMREX_SPACEDIM> force_slab_hi,
                                             Array<Real,AMREX_SPACEDIM> extforce,
                                             int update_massvel,int update_forces, amrex::Real mass_tolerance, int order_scheme)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();
    int extloads=external_loads_present;

    Real grav[]={AMREX_D_DECL(gravity[XDIR],gravity[YDIR],gravity[ZDIR])};
    Real slab_lo[]={AMREX_D_DECL(force_slab_lo[XDIR],force_slab_lo[YDIR],force_slab_lo[ZDIR])};
    Real slab_hi[]={AMREX_D_DECL(force_slab_hi[XDIR],force_slab_hi[YDIR],force_slab_hi[ZDIR])};
    Real extpforce[]={AMREX_D_DECL(extforce[XDIR],extforce[YDIR],extforce[ZDIR])};

    const int* lo = domain.loVect ();
    const int* hi = domain.hiVect ();


    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        //already nodal as mfi is from nodaldata
        const Box& nodalbox=mfi.validbox();

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        amrex::ParallelFor(nodalbox,[=]
        AMREX_GPU_DEVICE (int i,int j,int k) noexcept
        {
            if(update_massvel)
            {
                nodal_data_arr(i,j,k,MASS_INDEX)=zero;
                nodal_data_arr(i,j,k,VELX_INDEX)=zero;
                nodal_data_arr(i,j,k,VELY_INDEX)=zero;
                nodal_data_arr(i,j,k,VELZ_INDEX)=zero;
            }
            if(update_forces)
            {
                nodal_data_arr(i,j,k,FRCX_INDEX)=zero;
                nodal_data_arr(i,j,k,FRCY_INDEX)=zero;
                nodal_data_arr(i,j,k,FRCZ_INDEX)=zero;
            }
        });
    }

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


        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(nt,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
        	int lmin,lmax,nmin,nmax,mmin,mmax;

            ParticleType& p = pstruct[i];

            amrex::Real xp[AMREX_SPACEDIM];

            xp[XDIR]=p.pos(XDIR);
            xp[YDIR]=p.pos(YDIR);
            xp[ZDIR]=p.pos(ZDIR);

            auto iv = getParticleCell(p, plo, dxi, domain);

            if(order_scheme==1)
            {
            	lmin=0;
            	lmax=2;
            	nmin=0;
            	nmax=2;
            	mmin=0;
            	mmax=2;
            }
            else
            {
            	if(iv[0]==lo[0])
            	{
            		lmin=0;
            		lmax=3;
            	}
            	else if(iv[0]==hi[0])
            	{
            		lmin=-1;
            		lmax=2;
            	}
            	else
            	{
            		lmin=-1;
            		lmax=3;
            	}

            	if(iv[1]==lo[1])
            	{
            		mmin=0;
            		mmax=3;
            	}
            	else if(iv[1]==hi[1])
            	{
            		mmin=-1;
            	    mmax=2;
            	}
            	else
            	{
            		mmin=-1;
            		mmax=3;
            	}

            	if(iv[2]==lo[2])
            	{
            		nmin=0;
            		nmax=3;
            	}
            	else if(iv[2]==hi[2])
            	{
            		nmin=-1;
            		nmax=2;
            	}
            	else
            	{
            		nmin=-1;
            		nmax=3;
            	}
            }

            for(int n=nmin;n<nmax;n++)
            {
                for(int m=mmin;m<mmax;m++)
                {
                    for(int l=lmin;l<lmax;l++)
                    {
                        IntVect ivlocal(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n);

                        if(nodalbox.contains(ivlocal))
                        {

                            amrex::Real basisvalue=basisval(l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme,lo,hi);

                            if(update_massvel)
                            {

                                amrex::Real mass_contrib=p.rdata(realData::mass)*basisvalue;
                                amrex::Real p_contrib[AMREX_SPACEDIM] = 
                                {p.rdata(realData::mass)*p.rdata(realData::xvel)*basisvalue,
                                    p.rdata(realData::mass)*p.rdata(realData::yvel)*basisvalue,
                                    p.rdata(realData::mass)*p.rdata(realData::zvel)*basisvalue};


                                amrex::Gpu::Atomic::AddNoRet(&nodal_data_arr(ivlocal,MASS_INDEX), mass_contrib);

                                for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                                {
                                    amrex::Gpu::Atomic::AddNoRet(&nodal_data_arr(ivlocal,VELX_INDEX+dim),p_contrib[dim]);
                                }

                            }

                            if(update_forces)
                            {
                                amrex::Real basisval_grad[AMREX_SPACEDIM];
                                amrex::Real stress_tens[AMREX_SPACEDIM*AMREX_SPACEDIM];

                                get_tensor(p,realData::stress,stress_tens);

                                for(int d=0;d<AMREX_SPACEDIM;d++)
                                {
                                    basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,1,lo,hi);
                                }

                                amrex::Real bforce_contrib[AMREX_SPACEDIM]=
                                {   p.rdata(realData::mass)*grav[XDIR]*basisvalue,
                                    p.rdata(realData::mass)*grav[YDIR]*basisvalue,
                                    p.rdata(realData::mass)*grav[ZDIR]*basisvalue   };

                                if( extloads &&
                                    xp[XDIR]>force_slab_lo[XDIR] && xp[XDIR]<force_slab_hi[XDIR] &&
                                    xp[YDIR]>force_slab_lo[YDIR] && xp[YDIR]<force_slab_hi[YDIR] &&
                                    xp[ZDIR]>force_slab_lo[ZDIR] && xp[ZDIR]<force_slab_hi[ZDIR] )
                                {
                                    bforce_contrib[XDIR] += extpforce[XDIR]*basisvalue;
                                    bforce_contrib[YDIR] += extpforce[YDIR]*basisvalue;
                                    bforce_contrib[ZDIR] += extpforce[ZDIR]*basisvalue;
                                }

                                amrex::Real tensvect[AMREX_SPACEDIM];
                                tensor_vector_pdt(stress_tens,basisval_grad,tensvect);

                                amrex::Real intforce_contrib[AMREX_SPACEDIM]=
                                {-p.rdata(realData::volume)*tensvect[XDIR],
                                    -p.rdata(realData::volume)*tensvect[YDIR],
                                    -p.rdata(realData::volume)*tensvect[ZDIR]};

                                for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                                {
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,FRCX_INDEX+dim),
                                        bforce_contrib[dim]+intforce_contrib[dim]);
                                }
                            }
                        }

                    }
                }
            }
        });

    }
    //nodaldata.FillBoundary(geom.periodicity());
    //nodaldata.SumBoundary(geom.periodicity());
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

    	Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if(update_massvel)
                {
                    if(nodal_data_arr(i,j,k,MASS_INDEX) > 0.0)
                    {
                        for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                        {
                        	if(nodal_data_arr(i,j,k,MASS_INDEX)>=mass_tolerance)
                        	{
                        		nodal_data_arr(i,j,k,VELX_INDEX+dim)/=nodal_data_arr(i,j,k,MASS_INDEX);
                        	}
                        	else
                        	{
                        		nodal_data_arr(i,j,k,VELX_INDEX+dim) = 0.0;
                        	}
                        }
                    }
                }
            });

    }

}

void MPMParticleContainer::updatevolume(const amrex::Real& dt)
{
    BL_PROFILE("MPMParticleContainer::moveParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
    auto& plev  = GetParticles(lev);
    
    int periodic[AMREX_SPACEDIM]={Geom(lev).isPeriodic(XDIR),
                                  Geom(lev).isPeriodic(YDIR),
                                  Geom(lev).isPeriodic(ZDIR)};

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

            p.rdata(realData::jacobian) += (p.rdata(realData::strainrate+XX)+p.rdata(realData::strainrate+YY)+p.rdata(realData::strainrate+ZZ)) * dt * p.rdata(realData::jacobian);
            p.rdata(realData::volume)	= p.rdata(realData::vol_init)*p.rdata(realData::jacobian);
            p.rdata(realData::density)	= p.rdata(realData::mass)/p.rdata(realData::volume);

        });
    }
}

void MPMParticleContainer::moveParticles(const amrex::Real& dt)
{
    BL_PROFILE("MPMParticleContainer::moveParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
    auto& plev  = GetParticles(lev);

    int periodic[AMREX_SPACEDIM]={Geom(lev).isPeriodic(XDIR),
                                  Geom(lev).isPeriodic(YDIR),
                                  Geom(lev).isPeriodic(ZDIR)};

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

            p.pos(0) += p.rdata(realData::xvel) * dt;
            p.pos(1) += p.rdata(realData::yvel) * dt;
            p.pos(2) += p.rdata(realData::zvel) * dt;

            if (!periodic[XDIR] && (p.pos(0) < plo[0]))
            {
                p.pos(0) = two*plo[0] - p.pos(0);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (!periodic[XDIR] && (p.pos(0) > phi[0]))
            {
                p.pos(0) = two*phi[0] - p.pos(0);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (!periodic[YDIR] && (p.pos(1) < plo[1]))
            {
                p.pos(1) = two*plo[1] - p.pos(1);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (!periodic[YDIR] && (p.pos(1) > phi[1]))
            {
                p.pos(1) = two*phi[1] - p.pos(1);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (!periodic[ZDIR] && (p.pos(2) < plo[2]))
            {
                p.pos(2) = two*plo[2] - p.pos(2);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }
            if (!periodic[ZDIR] && (p.pos(2) > phi[2]))
            {
                p.pos(2) = two*phi[2] - p.pos(2);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }

        });
    }
}


void MPMParticleContainer::interpolate_mass_from_grid(MultiFab& nodaldata,int order_scheme)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    int ncomp=nodaldata.nComp();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const int np = aos.numRealParticles();

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];

            amrex::Real xp[AMREX_SPACEDIM];
            amrex::Real gradvp[AMREX_SPACEDIM][AMREX_SPACEDIM]={0.0};

            xp[XDIR]=p.pos(XDIR);
            xp[YDIR]=p.pos(YDIR);
            xp[ZDIR]=p.pos(ZDIR);

            auto iv = getParticleCell(p, plo, dxi, domain);

            if(order_scheme==1)
            {
            	p.rdata(realData::vol_init) = bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,MASS_INDEX);
            	p.rdata(realData::vol_init) = p.rdata(realData::vol_init)/(dx[XDIR]*dx[YDIR]*dx[ZDIR]);//Actually the density. Dont get confused by the variable name.
            	p.rdata(realData::vol_init) = p.rdata(realData::mass)/p.rdata(realData::vol_init);	//INitial volume
            }
        });
    }
}


void MPMParticleContainer::move_particles_from_nodevel(MultiFab& nodaldata,const amrex::Real& dt,int order_scheme)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto domain = geom.Domain();

    int periodic[AMREX_SPACEDIM]={Geom(lev).isPeriodic(XDIR),
                                      Geom(lev).isPeriodic(YDIR),
                                      Geom(lev).isPeriodic(ZDIR)};

    int ncomp=nodaldata.nComp();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const int np = aos.numRealParticles();

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];

            amrex::Real xp[AMREX_SPACEDIM];

            xp[XDIR]=p.pos(XDIR);
            xp[YDIR]=p.pos(YDIR);
            xp[ZDIR]=p.pos(ZDIR);

            auto iv = getParticleCell(p, plo, dxi, domain);
           	if(order_scheme==1)
           	{
           		p.pos(XDIR) += bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELX_INDEX)*dt;
           		p.pos(YDIR) += bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELY_INDEX)*dt;
           		p.pos(ZDIR) += bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELZ_INDEX)*dt;
            }

           	if (!periodic[XDIR] && (p.pos(0) < plo[0]))
           	{
           		p.pos(0) = two*plo[0] - p.pos(0);
           		p.rdata(realData::xvel) = -p.rdata(realData::xvel);
           	}
           	if (!periodic[XDIR] && (p.pos(0) > phi[0]))
           	{
           		p.pos(0) = two*phi[0] - p.pos(0);
           		p.rdata(realData::xvel) = -p.rdata(realData::xvel);
           	}
           	if (!periodic[YDIR] && (p.pos(1) < plo[1]))
           	{
           		p.pos(1) = two*plo[1] - p.pos(1);
           		p.rdata(realData::yvel) = -p.rdata(realData::yvel);
           	}
           	if (!periodic[YDIR] && (p.pos(1) > phi[1]))
           	{
           		p.pos(1) = two*phi[1] - p.pos(1);
           		p.rdata(realData::yvel) = -p.rdata(realData::yvel);
           	}
           	if (!periodic[ZDIR] && (p.pos(2) < plo[2]))
           	{
           		p.pos(2) = two*plo[2] - p.pos(2);
           		p.rdata(realData::zvel) = -p.rdata(realData::zvel);
           	}
           	if (!periodic[ZDIR] && (p.pos(2) > phi[2]))
           	{
           		p.pos(2) = two*phi[2] - p.pos(2);
           		p.rdata(realData::zvel) = -p.rdata(realData::zvel);
           	}
        });
    }
}



void MPMParticleContainer::interpolate_from_grid(MultiFab& nodaldata,int update_vel,int update_strainrate,int order_scheme,amrex::Real alpha_pic_flip)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    int ncomp=nodaldata.nComp();
    const int* lo = domain.loVect ();
    const int* hi = domain.hiVect ();

    nodaldata.FillBoundary(geom.periodicity());

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const int np = aos.numRealParticles();

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
        	int lmin,lmax,nmin,nmax,mmin,mmax;
            ParticleType& p = pstruct[i];

            amrex::Real xp[AMREX_SPACEDIM];
            amrex::Real gradvp[AMREX_SPACEDIM][AMREX_SPACEDIM]={0.0};

            xp[XDIR]=p.pos(XDIR);
            xp[YDIR]=p.pos(YDIR);
            xp[ZDIR]=p.pos(ZDIR);

            auto iv = getParticleCell(p, plo, dxi, domain);

            if(order_scheme==1)
            {
            	lmin=0;
            	lmax=2;
            	nmin=0;
            	nmax=2;
            	mmin=0;
            	mmax=2;
            }
            else if(order_scheme==3)
            {
            	if(iv[0]==lo[0])
            	{
            		lmin=0;
            		lmax=3;
            	}
            	else if(iv[0]==hi[0])
            	{
            		lmin=-1;
            		lmax=2;
            	}
            	else
            	{
            		lmin=-1;
            		lmax=3;
            	}
            	if(iv[1]==lo[1])
            	{
            		mmin=0;
            		mmax=3;
            	}
            	else if(iv[1]==hi[1])
            	{
            		mmin=-1;
            		mmax=2;
            	}
            	else
            	{
            		mmin=-1;
            		mmax=3;
            	}
            	if(iv[2]==lo[2])
            	{
            		nmin=0;
            		nmax=3;
            	}
            	else if(iv[2]==hi[2])
            	{
            		nmin=-1;
            		nmax=2;
            	}
            	else
            	{
            		nmin=-1;
            		nmax=3;
            	}
            }


            if(update_vel)
            {
            	if(order_scheme==1)
            	{
            		p.rdata(realData::xvel) = (alpha_pic_flip)*p.rdata(realData::xvel)
            								  +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELX_INDEX)
											  +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELX_INDEX);
            		p.rdata(realData::yvel) = (alpha_pic_flip)*p.rdata(realData::yvel)
            								  +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELY_INDEX)
											  +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELY_INDEX);
            		p.rdata(realData::zvel) = (alpha_pic_flip)*p.rdata(realData::zvel)
            								  +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELZ_INDEX)
											  +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELZ_INDEX);
            	}
            	if(order_scheme==3)
            	{
            		/*p.rdata(realData::xvel) = (alpha_pic_flip)*p.rdata(realData::xvel)
            								  +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELX_INDEX)
											  +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELX_INDEX);
            		p.rdata(realData::yvel) = (alpha_pic_flip)*p.rdata(realData::yvel)
            								  +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELY_INDEX)
											  +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELY_INDEX);
            		p.rdata(realData::zvel) = (alpha_pic_flip)*p.rdata(realData::zvel)
            								  +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELZ_INDEX)
											  +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELZ_INDEX);*/


            		p.rdata(realData::xvel) = (alpha_pic_flip)*p.rdata(realData::xvel)
											  +(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELX_INDEX,lo,hi)
											  +(1-alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELX_INDEX,lo,hi);
            		p.rdata(realData::yvel) = (alpha_pic_flip)*p.rdata(realData::yvel)
											  +(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELY_INDEX,lo,hi)
											  +(1-alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELY_INDEX,lo,hi);
            		p.rdata(realData::zvel) = (alpha_pic_flip)*p.rdata(realData::zvel)
											  +(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELZ_INDEX,lo,hi)
											  +(1-alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELZ_INDEX,lo,hi);


            		/*
            		p.rdata(realData::xvel) = cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELX_INDEX,lo,hi);
            		p.rdata(realData::yvel) = cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELY_INDEX,lo,hi);
            		p.rdata(realData::zvel) = cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELZ_INDEX,lo,hi);
            		*/

            	}
            }

            if(update_strainrate)
            {
                for(int n=nmin;n<nmax;n++)
                {
                    for(int m=mmin;m<mmax;m++)
                    {
                        for(int l=lmin;l<lmax;l++)
                        {
                            amrex::Real basisval_grad[AMREX_SPACEDIM];
                            for(int d=0;d<AMREX_SPACEDIM;d++)
                            {
                                basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme,lo,hi);
                            }

                            gradvp[XDIR][XDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELX_INDEX)*basisval_grad[XDIR];
                            gradvp[XDIR][YDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELX_INDEX)*basisval_grad[YDIR];
                            gradvp[XDIR][ZDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELX_INDEX)*basisval_grad[ZDIR];

                            gradvp[YDIR][XDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELY_INDEX)*basisval_grad[XDIR];
                            gradvp[YDIR][YDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELY_INDEX)*basisval_grad[YDIR];
                            gradvp[YDIR][ZDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELY_INDEX)*basisval_grad[ZDIR];

                            gradvp[ZDIR][XDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELZ_INDEX)*basisval_grad[XDIR];
                            gradvp[ZDIR][YDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELZ_INDEX)*basisval_grad[YDIR];
                            gradvp[ZDIR][ZDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELZ_INDEX)*basisval_grad[ZDIR];

                        }
                    }
                }

                int ind=0;
                for(int d1=0;d1<AMREX_SPACEDIM;d1++)
                {
                    for(int d2=d1;d2<AMREX_SPACEDIM;d2++)
                    {
                        p.rdata(realData::strainrate+ind)=0.5*(gradvp[d1][d2]+gradvp[d2][d1]);
                        ind++;
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

    WritePlotFile(pltfile, "particles",writeflags_real, 
                  writeflags_int, real_data_names, int_data_names);
}
