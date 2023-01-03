#include <mpm_particle_container.H>
#include <interpolants.H>

int MPMParticleContainer::checkifrigidnodespresent()
{
	int rigidnodespresent=0;
	const int lev = 0;
	auto& plev  = GetParticles(lev);

	using PType = typename MPMParticleContainer::SuperParticleType;
	rigidnodespresent = amrex::ReduceMax(*this, [=]
	    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
	    {
	        if(p.idata(intData::phase)==1)
	        {
	        	int rigidnodespresenttmp=1;
				return(rigidnodespresenttmp);
	        }
	        else
	        {
	        	int rigidnodespresenttmp=0;
	        	return(rigidnodespresenttmp);
	        }
	    });

	#ifdef BL_USE_MPI
	    ParallelDescriptor::ReduceIntMax(rigidnodespresent);
	#endif
	return(rigidnodespresent);

}

int MPMParticleContainer::get_num_of_bodies()
{
	int rigidnodespresent=0;
	const int lev = 0;
	auto& plev  = GetParticles(lev);

	int num_of_body=0;
	using PType = typename MPMParticleContainer::SuperParticleType;
	num_of_body = amrex::ReduceMax(*this, [=]
	    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
	    {
	        if(p.idata(intData::phase)==1)
	        {
	        	int rigidnodespresenttmp=p.idata(intData::rigid_body_id);
				return(rigidnodespresenttmp);
	        }
	        else
	        {
	        	int rigidnodespresenttmp=0;
	        	return(rigidnodespresenttmp);
	        }
	    });

	#ifdef BL_USE_MPI
	    ParallelDescriptor::ReduceIntMax(num_of_body);
	#endif
	return(num_of_body);

}

int MPMParticleContainer::Calculate_Total_Number_of_rigid_particles(int body_id)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    int total_num=0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_num = amrex::ReduceSum(*this, [=]
        AMREX_GPU_HOST_DEVICE (const PType& p) -> int
        {
    		if(p.idata(intData::phase)==1 and p.idata(intData::rigid_body_id)==body_id)
    		{
    			int num=1;
    			return(num);
    		}
    		else
    		{
    			int num=0;
    			return(num);
    		}
        });

	#ifdef BL_USE_MPI
	    ParallelDescriptor::ReduceIntSum(total_num);
	#endif

	    return(total_num);
}

void MPMParticleContainer::Calculate_Total_Number_of_MaterialParticles(int &total_num)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    total_num=0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_num = amrex::ReduceSum(*this, [=]
        AMREX_GPU_HOST_DEVICE (const PType& p) -> int
        {
    		if(p.idata(intData::phase)==0)
    		{
    			return(1);
    		}
    		else
    		{
    			return(0);
    		}
        });

	#ifdef BL_USE_MPI
	    ParallelDescriptor::ReduceIntSum(total_num);
	#endif
}

void MPMParticleContainer::Calculate_Total_Mass_RigidParticles(int body_id,Real &total_mass)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    total_mass=0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_mass = amrex::ReduceSum(*this, [=]
        AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
    		if(p.idata(intData::phase)==1 and p.idata(intData::rigid_body_id)==body_id)
    		{
    			return(p.rdata(realData::mass));
    		}
    		else
    		{
    			return(0.0);
    		}
        });

	#ifdef BL_USE_MPI
	    ParallelDescriptor::ReduceRealSum(total_mass);
	#endif
}

void MPMParticleContainer::Calculate_Total_Mass_MaterialPoints(Real &total_mass)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    total_mass=0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_mass = amrex::ReduceSum(*this, [=]
        AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
    		if(p.idata(intData::phase)==0)
    		{
    			return(p.rdata(realData::mass));
    		}
    		else
    		{
    			return(0.0);
    		}
        });

	#ifdef BL_USE_MPI
	    ParallelDescriptor::ReduceRealSum(total_mass);
	#endif
}

void MPMParticleContainer::Calculate_Total_Vol_MaterialPoints(Real &total_vol)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    total_vol=0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_vol = amrex::ReduceSum(*this, [=]
        AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
    		if(p.idata(intData::phase)==0)
    		{
    			return(p.rdata(realData::volume));
    		}
    		else
    		{
    			return(0.0);
    		}
        });

	#ifdef BL_USE_MPI
	    ParallelDescriptor::ReduceRealSum(total_vol);
	#endif
}

amrex::Real MPMParticleContainer::Calculate_Total_Vol_RigidParticles(int body_id)
{

	amrex::Real total_vol=0.0;
	const int lev = 0;
	auto& plev  = GetParticles(lev);

	using PType = typename MPMParticleContainer::SuperParticleType;
	total_vol = amrex::ReduceSum(*this, [=]
	    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
	    {
			if(p.idata(intData::phase)==1 and  p.idata(intData::rigid_body_id)==body_id)
	        {
	        	return(p.rdata(realData::volume));
	        }
			else
			{
				return(0.0);
			}
	    });

	#ifdef BL_USE_MPI
	    ParallelDescriptor::ReduceRealSum(total_vol);
	#endif
	return(total_vol);
}

void MPMParticleContainer::deposit_onto_grid(MultiFab& nodaldata,
                                             Array<Real,AMREX_SPACEDIM> gravity,
                                             int external_loads_present,
                                             Array<Real,AMREX_SPACEDIM> force_slab_lo,
                                             Array<Real,AMREX_SPACEDIM> force_slab_hi,
                                             Array<Real,AMREX_SPACEDIM> extforce,
                                             int update_massvel,
				             				 int update_forces,
					 						 amrex::Real mass_tolerance,
					     					 GpuArray<int,AMREX_SPACEDIM> order_scheme_directional,
											 GpuArray<int,AMREX_SPACEDIM> periodic)
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
    
    const int* loarr = domain.loVect ();
    const int* hiarr = domain.hiVect ();
    
    int lo[]={loarr[0],loarr[1],loarr[2]};
    int hi[]={hiarr[0],hiarr[1],hiarr[2]};

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
            if(update_forces==2)
            {
            	nodal_data_arr(i,j,k,STRESS_INDEX)=zero;
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

            if(p.idata(intData::phase)==0)		//Compute only for standard particles and not rigid particles with phase=1
            {

            	amrex::Real xp[AMREX_SPACEDIM];

            	xp[XDIR]=p.pos(XDIR);
            	xp[YDIR]=p.pos(YDIR);
            	xp[ZDIR]=p.pos(ZDIR);

            	auto iv = getParticleCell(p, plo, dxi, domain);

            	lmin=(order_scheme_directional[0]==1)?0:((order_scheme_directional[0]==3)?(iv[XDIR]==lo[XDIR])?0:((iv[XDIR]==hi[XDIR])?-1:-1):-1000);
            	lmax=(order_scheme_directional[0]==1)?2:((order_scheme_directional[0]==3)?(iv[XDIR]==lo[XDIR])?lmin+3:((iv[XDIR]==hi[XDIR])?lmin+3:lmin+4):-1000);

            	mmin=(order_scheme_directional[1]==1)?0:((order_scheme_directional[1]==3)?(iv[YDIR]==lo[YDIR])?0:((iv[YDIR]==hi[YDIR])?-1:-1):-1000);
            	mmax=(order_scheme_directional[1]==1)?2:((order_scheme_directional[1]==3)?(iv[YDIR]==lo[YDIR])?mmin+3:((iv[YDIR]==hi[YDIR])?mmin+3:mmin+4):-1000);

            	nmin=(order_scheme_directional[2]==1)?0:((order_scheme_directional[2]==3)?(iv[ZDIR]==lo[ZDIR])?0:((iv[ZDIR]==hi[ZDIR])?-1:-1):-1000);
            	nmax=(order_scheme_directional[2]==1)?2:((order_scheme_directional[2]==3)?(iv[ZDIR]==lo[ZDIR])?nmin+3:((iv[ZDIR]==hi[ZDIR])?nmin+3:nmin+4):-1000);


            	if(lmin==-1000 or lmax==-1000 or mmin==-1000 or mmax==-1000 or nmin==-1000 or nmax==-1000)
            	{
            		amrex::Abort("\nError. Something wrong with min/max index values in deposit onto grid");
            	}

            	for(int n=nmin;n<nmax;n++)
            	{
            		for(int m=mmin;m<mmax;m++)
            		{
            			for(int l=lmin;l<lmax;l++)
            			{
            				IntVect ivlocal(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n);

            				if(iv[YDIR]+m==lo[1] && iv[XDIR]+l==25 && iv[ZDIR]+n==25)
            				//if(iv[XDIR]+l==25 && iv[ZDIR]+n==25)
            				{
            					//amrex::Print()<<"\n Particle = "<<p.pos(0)<<" "<<p.pos(1)<<" "<<p.pos(2);
            				}
            				if(nodalbox.contains(ivlocal))
            				{
            					amrex::Real basisvalue=basisval(l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme_directional,periodic,lo,hi);

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
            							basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme_directional,periodic,lo,hi);
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
            						{	-p.rdata(realData::volume)*tensvect[XDIR],
            							-p.rdata(realData::volume)*tensvect[YDIR],
										-p.rdata(realData::volume)*tensvect[ZDIR]};

            						for(int dim=0;dim<AMREX_SPACEDIM;dim++)
            						{
            							amrex::Gpu::Atomic::AddNoRet(
            									&nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,FRCX_INDEX+dim),
												bforce_contrib[dim]+intforce_contrib[dim]);
            						}
            					}

            					if(update_forces==2)
            					{
            						amrex::Real stress_contrib=p.rdata(realData::stress+3)*p.rdata(realData::mass)*basisvalue;
            						amrex::Gpu::Atomic::AddNoRet(&nodal_data_arr(ivlocal,STRESS_INDEX), stress_contrib);
            					}
            				} //nodalbox if loop
            			} //l loop
            		}//m loop
            	}//n loop
            }//phase=0 if loop
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
            	//amrex::Print()<<"\n Nodal mass values for i = "<<i<<" j = "<<j<<" k = "<<k<<" is "<<nodal_data_arr(i,j,k,MASS_INDEX);
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
            if(update_forces==2)
            {
            	if(nodal_data_arr(i,j,k,MASS_INDEX) > 0.0)
            	{
            		if(nodal_data_arr(i,j,k,MASS_INDEX)>=mass_tolerance)
            		{
            			nodal_data_arr(i,j,k,STRESS_INDEX)/=nodal_data_arr(i,j,k,MASS_INDEX);
            		}
            		else
            		{
            			nodal_data_arr(i,j,k,STRESS_INDEX) = 0.0;
            		}
            	}
            }
        });

    }

}




void MPMParticleContainer::deposit_onto_grid_rigidnodesonly(MultiFab& nodaldata,
                                             Array<Real,AMREX_SPACEDIM> gravity,
                                             int external_loads_present,
                                             Array<Real,AMREX_SPACEDIM> force_slab_lo,
                                             Array<Real,AMREX_SPACEDIM> force_slab_hi,
                                             Array<Real,AMREX_SPACEDIM> extforce,
                                             int update_massvel,int update_forces, amrex::Real mass_tolerance,
                                             GpuArray<int,AMREX_SPACEDIM> order_scheme_directional,
					     GpuArray<int,AMREX_SPACEDIM> periodic)
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

    const int* loarr = domain.loVect ();
    const int* hiarr = domain.hiVect ();
    
    int lo[]={loarr[0],loarr[1],loarr[2]};
    int hi[]={hiarr[0],hiarr[1],hiarr[2]};

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
    	//already nodal as mfi is from nodaldata
    	const Box& nodalbox=mfi.validbox();

    	Array4<Real> nodal_data_arr=nodaldata.array(mfi);

    	amrex::ParallelFor(nodalbox,[=]
			AMREX_GPU_DEVICE (int i,int j,int k) noexcept
            {
            	nodal_data_arr(i,j,k,RIGID_BODY_ID)=-1;
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


            if(p.idata(intData::phase)==1)		//Compute only for rigid particles with phase=1
            {
            	amrex::Real xp[AMREX_SPACEDIM];

            	xp[XDIR]=p.pos(XDIR);
            	xp[YDIR]=p.pos(YDIR);
            	xp[ZDIR]=p.pos(ZDIR);

            	auto iv = getParticleCell(p, plo, dxi, domain);

            	lmin=(order_scheme_directional[0]==1)?0:((order_scheme_directional[0]==3)?(iv[XDIR]==lo[XDIR])?0:((iv[XDIR]==hi[XDIR])?-1:-1):-1000);
            	lmax=(order_scheme_directional[0]==1)?2:((order_scheme_directional[0]==3)?(iv[XDIR]==lo[XDIR])?lmin+3:((iv[XDIR]==hi[XDIR])?lmin+3:lmin+4):-1000);

            	mmin=(order_scheme_directional[1]==1)?0:((order_scheme_directional[1]==3)?(iv[YDIR]==lo[YDIR])?0:((iv[YDIR]==hi[YDIR])?-1:-1):-1000);
            	mmax=(order_scheme_directional[1]==1)?2:((order_scheme_directional[1]==3)?(iv[YDIR]==lo[YDIR])?mmin+3:((iv[YDIR]==hi[YDIR])?mmin+3:mmin+4):-1000);

            	nmin=(order_scheme_directional[2]==1)?0:((order_scheme_directional[2]==3)?(iv[ZDIR]==lo[ZDIR])?0:((iv[ZDIR]==hi[ZDIR])?-1:-1):-1000);
            	nmax=(order_scheme_directional[2]==1)?2:((order_scheme_directional[2]==3)?(iv[ZDIR]==lo[ZDIR])?nmin+3:((iv[ZDIR]==hi[ZDIR])?nmin+3:nmin+4):-1000);


            	if(lmin==-1000 or lmax==-1000 or mmin==-1000 or mmax==-1000 or nmin==-1000 or nmax==-1000)
            	{
            		amrex::Abort("\nError. Something wrong with min/max index values in deposit onto grid");
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

            					amrex::Real basisvalue=basisval(l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme_directional,periodic,lo,hi);

            						amrex::Real mass_contrib=p.rdata(realData::mass)*basisvalue;
            						amrex::Real p_contrib[AMREX_SPACEDIM] =
            						{p.rdata(realData::mass)*p.rdata(realData::xvel)*basisvalue,
                                     p.rdata(realData::mass)*p.rdata(realData::yvel)*basisvalue,
                                     p.rdata(realData::mass)*p.rdata(realData::zvel)*basisvalue};

            						amrex::Gpu::Atomic::AddNoRet(&nodal_data_arr(ivlocal,MASS_RIGID_INDEX), mass_contrib);
            						nodal_data_arr(ivlocal,RIGID_BODY_ID)=p.idata(intData::rigid_body_id);

            						for(int dim=0;dim<AMREX_SPACEDIM;dim++)
            						{
            							amrex::Gpu::Atomic::AddNoRet(&nodal_data_arr(ivlocal,VELX_RIGID_INDEX+dim),p_contrib[dim]);
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
            	//amrex::Print()<<"\n Nodal mass values for i = "<<i<<" j = "<<j<<" k = "<<k<<" is "<<nodal_data_arr(i,j,k,MASS_INDEX);
                if(nodal_data_arr(i,j,k,MASS_RIGID_INDEX) > 0.0)
                {
                    for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                    {
                        if(nodal_data_arr(i,j,k,MASS_RIGID_INDEX)>=mass_tolerance)
                        {
                            nodal_data_arr(i,j,k,VELX_RIGID_INDEX+dim)/=nodal_data_arr(i,j,k,MASS_RIGID_INDEX);
                        }
                        else
                        {
                            nodal_data_arr(i,j,k,VELX_RIGID_INDEX+dim) = 0.0;
                        }
                    }
                }
            }
        });

    }

}

void MPMParticleContainer::interpolate_from_grid(MultiFab& nodaldata,int update_vel,
                    int update_strainrate,
					GpuArray <int,AMREX_SPACEDIM> order_scheme_directional,
					GpuArray <int,AMREX_SPACEDIM> periodic,
					amrex::Real alpha_pic_flip,
					amrex::Real dt)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    int ncomp=nodaldata.nComp();
    const int* loarr = domain.loVect ();
    const int* hiarr = domain.hiVect ();
    
    int lo[]={loarr[0],loarr[1],loarr[2]};
    int hi[]={hiarr[0],hiarr[1],hiarr[2]};
    const double pi = 3.141592654;

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

            if(p.idata(intData::phase)==0)
            {

				amrex::Real xp[AMREX_SPACEDIM];
				amrex::Real gradvp[AMREX_SPACEDIM][AMREX_SPACEDIM]={0.0};

				xp[XDIR]=p.pos(XDIR);
				xp[YDIR]=p.pos(YDIR);
				xp[ZDIR]=p.pos(ZDIR);

				auto iv = getParticleCell(p, plo, dxi, domain);

				lmin=(order_scheme_directional[0]==1)?0:((order_scheme_directional[0]==3)?(iv[XDIR]==lo[XDIR])?0:((iv[XDIR]==hi[XDIR])?-1:-1):-1000);
				lmax=(order_scheme_directional[0]==1)?2:((order_scheme_directional[0]==3)?(iv[XDIR]==lo[XDIR])?lmin+3:((iv[XDIR]==hi[XDIR])?lmin+3:lmin+4):-1000);

				mmin=(order_scheme_directional[1]==1)?0:((order_scheme_directional[1]==3)?(iv[YDIR]==lo[YDIR])?0:((iv[YDIR]==hi[YDIR])?-1:-1):-1000);
				mmax=(order_scheme_directional[1]==1)?2:((order_scheme_directional[1]==3)?(iv[YDIR]==lo[YDIR])?mmin+3:((iv[YDIR]==hi[YDIR])?mmin+3:mmin+4):-1000);

				nmin=(order_scheme_directional[2]==1)?0:((order_scheme_directional[2]==3)?(iv[ZDIR]==lo[ZDIR])?0:((iv[ZDIR]==hi[ZDIR])?-1:-1):-1000);
				nmax=(order_scheme_directional[2]==1)?2:((order_scheme_directional[2]==3)?(iv[ZDIR]==lo[ZDIR])?nmin+3:((iv[ZDIR]==hi[ZDIR])?nmin+3:nmin+4):-1000);

				if(lmin ==-1000 or lmax==-1000 or mmin==-1000 or mmax==-1000 or nmin==-1000 or nmax==-1000)
				{
					amrex::Abort("\nError. Something wrong with min/max index values");
				}

				if(update_vel)
				{
					if(order_scheme_directional[0]==1)
					{

						p.rdata(realData::xvel_prime) = bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELX_INDEX);
						p.rdata(realData::xvel) = (alpha_pic_flip)*p.rdata(realData::xvel)
						+(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELX_INDEX)
						+(1-alpha_pic_flip)*p.rdata(realData::xvel_prime);
					}
					else if(order_scheme_directional[0]==3)
					{
						p.rdata(realData::xvel_prime) = cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELX_INDEX,lo,hi);
						p.rdata(realData::xvel) = (alpha_pic_flip)*p.rdata(realData::xvel)
						+(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELX_INDEX,lo,hi)
						+(1-alpha_pic_flip)*p.rdata(realData::xvel_prime);
					}

					if(order_scheme_directional[1]==1)
					{
						p.rdata(realData::yacceleration)= p.rdata(realData::yvel);
						p.rdata(realData::yvel_prime) = bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELY_INDEX);
						p.rdata(realData::yvel) = (alpha_pic_flip)*p.rdata(realData::yvel)
						+(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELY_INDEX)
						+(1-alpha_pic_flip)*p.rdata(realData::yvel_prime);
						p.rdata(realData::yacceleration)= (p.rdata(realData::yvel)-p.rdata(realData::yacceleration))/dt;
					}
					else if(order_scheme_directional[1]==3)
					{
						p.rdata(realData::yacceleration)= p.rdata(realData::yvel);
						p.rdata(realData::yvel_prime) = cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELY_INDEX,lo,hi);
						p.rdata(realData::yvel) = (alpha_pic_flip)*p.rdata(realData::yvel)
						+(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELY_INDEX,lo,hi)
						+(1-alpha_pic_flip)*p.rdata(realData::yvel_prime);
						p.rdata(realData::yacceleration)= (p.rdata(realData::yvel)-p.rdata(realData::yacceleration))/dt;
					}

					if(order_scheme_directional[2]==1)
					{
						p.rdata(realData::zvel_prime) = bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELZ_INDEX);
						p.rdata(realData::zvel) = (alpha_pic_flip)*p.rdata(realData::zvel)
						+(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELZ_INDEX)
						+(1-alpha_pic_flip)*p.rdata(realData::zvel_prime);
					}
					else if(order_scheme_directional[2]==3)
					{
						p.rdata(realData::zvel_prime) = cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELZ_INDEX,lo,hi);
						p.rdata(realData::zvel) = (alpha_pic_flip)*p.rdata(realData::zvel)
						+(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELZ_INDEX,lo,hi)
						+(1-alpha_pic_flip)*p.rdata(realData::zvel_prime);
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
									basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme_directional,periodic,lo,hi);
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

					//Calculate deformation gradient tensor at time t+dt
					get_deformation_gradient_tensor(p,realData::deformation_gradient,gradvp,dt);

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
            }
        });
    }
}


void MPMParticleContainer::calculate_nodal_normal	(MultiFab& nodaldata,
														 amrex::Real mass_tolerance,
														 GpuArray<int,AMREX_SPACEDIM> order_scheme_directional,
														 GpuArray<int,AMREX_SPACEDIM> periodic)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    const int* loarr = domain.loVect ();
    const int* hiarr = domain.hiVect ();

    int lo[]={loarr[0],loarr[1],loarr[2]};
    int hi[]={hiarr[0],hiarr[1],hiarr[2]};

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        const Box& nodalbox=mfi.validbox();

        Array4<Real> nodal_data_arr=nodaldata.array(mfi);

        amrex::ParallelFor(nodalbox,[=]
        AMREX_GPU_DEVICE (int i,int j,int k) noexcept
        {
        	nodal_data_arr(i,j,k,NORMALX)=zero;
        	nodal_data_arr(i,j,k,NORMALY)=zero;
        	nodal_data_arr(i,j,k,NORMALZ)=zero;
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

            if(p.idata(intData::phase)==0)		//Compute only for standard particles and not rigid particles with phase=1
            {

            	amrex::Real xp[AMREX_SPACEDIM];

            	xp[XDIR]=p.pos(XDIR);
            	xp[YDIR]=p.pos(YDIR);
            	xp[ZDIR]=p.pos(ZDIR);

            	auto iv = getParticleCell(p, plo, dxi, domain);

            	lmin=(order_scheme_directional[0]==1)?0:((order_scheme_directional[0]==3)?(iv[XDIR]==lo[XDIR])?0:((iv[XDIR]==hi[XDIR])?-1:-1):-1000);
            	lmax=(order_scheme_directional[0]==1)?2:((order_scheme_directional[0]==3)?(iv[XDIR]==lo[XDIR])?lmin+3:((iv[XDIR]==hi[XDIR])?lmin+3:lmin+4):-1000);

            	mmin=(order_scheme_directional[1]==1)?0:((order_scheme_directional[1]==3)?(iv[YDIR]==lo[YDIR])?0:((iv[YDIR]==hi[YDIR])?-1:-1):-1000);
            	mmax=(order_scheme_directional[1]==1)?2:((order_scheme_directional[1]==3)?(iv[YDIR]==lo[YDIR])?mmin+3:((iv[YDIR]==hi[YDIR])?mmin+3:mmin+4):-1000);

            	nmin=(order_scheme_directional[2]==1)?0:((order_scheme_directional[2]==3)?(iv[ZDIR]==lo[ZDIR])?0:((iv[ZDIR]==hi[ZDIR])?-1:-1):-1000);
            	nmax=(order_scheme_directional[2]==1)?2:((order_scheme_directional[2]==3)?(iv[ZDIR]==lo[ZDIR])?nmin+3:((iv[ZDIR]==hi[ZDIR])?nmin+3:nmin+4):-1000);


            	if(lmin==-1000 or lmax==-1000 or mmin==-1000 or mmax==-1000 or nmin==-1000 or nmax==-1000)
            	{
            		amrex::Abort("\nError. Something wrong with min/max index values in deposit onto grid");
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

            					amrex::Real basisval_grad[AMREX_SPACEDIM];
            					for(int d=0;d<AMREX_SPACEDIM;d++)
            					{
            						basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme_directional,periodic,lo,hi);
            					}
            					amrex::Real normal[AMREX_SPACEDIM]={p.rdata(realData::mass)*basisval_grad[XDIR],p.rdata(realData::mass)*basisval_grad[YDIR],p.rdata(realData::mass)*basisval_grad[ZDIR]};
            					for(int dim=0;dim<AMREX_SPACEDIM;dim++)
            					{
            						amrex::Gpu::Atomic::AddNoRet(&nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,NORMALX+dim),normal[dim]);
            					}
            				}
            			}
            		}
            	}
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

        amrex::ParallelFor(
        nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
        	amrex::Real nmag=pow((nodal_data_arr(i,j,k,NORMALX)*nodal_data_arr(i,j,k,NORMALX)+nodal_data_arr(i,j,k,NORMALY)*nodal_data_arr(i,j,k,NORMALY)+nodal_data_arr(i,j,k,NORMALZ)*nodal_data_arr(i,j,k,NORMALZ)),0.5);

        	if(nmag>mass_tolerance)
        	{
        		for(int d=0;d<AMREX_SPACEDIM;d++)
        		{
        			nodal_data_arr(i,j,k,NORMALX+d)=nodal_data_arr(i,j,k,NORMALX+d)/nmag;
        		}
        	}
        	else
        	{
        		nodal_data_arr(i,j,k,NORMALX)=0.0;
        		nodal_data_arr(i,j,k,NORMALY)=0.0;
        		nodal_data_arr(i,j,k,NORMALZ)=0.0;

        	}
        });
    }
}
