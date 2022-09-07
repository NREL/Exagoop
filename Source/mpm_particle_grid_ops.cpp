#include <mpm_particle_container.H>
#include <interpolants.H>

void MPMParticleContainer::deposit_onto_grid(MultiFab& nodaldata,
                                             Array<Real,AMREX_SPACEDIM> gravity,
                                             int external_loads_present,
                                             Array<Real,AMREX_SPACEDIM> force_slab_lo,
                                             Array<Real,AMREX_SPACEDIM> force_slab_hi,
                                             Array<Real,AMREX_SPACEDIM> extforce,
                                             int update_massvel,int update_forces, amrex::Real mass_tolerance,
											 Array<int,AMREX_SPACEDIM> order_scheme_directional,
											 Array<int,AMREX_SPACEDIM> periodic)
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

            p.rdata(realData::vol_init) = bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,MASS_INDEX);
                
            //Actually the density. Dont get confused by the variable name.
            p.rdata(realData::vol_init) = p.rdata(realData::vol_init)/(dx[XDIR]*dx[YDIR]*dx[ZDIR]);
                
            //INitial volume
            p.rdata(realData::vol_init) = p.rdata(realData::mass)/p.rdata(realData::vol_init);

        });
    }
}

void MPMParticleContainer::interpolate_from_grid(MultiFab& nodaldata,int update_vel,
                    int update_strainrate,
					Array <int,AMREX_SPACEDIM> order_scheme_directional,
					Array <int,AMREX_SPACEDIM> periodic,
					amrex::Real alpha_pic_flip)
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
                    p.rdata(realData::xvel) = (alpha_pic_flip)*p.rdata(realData::xvel)
                    +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELX_INDEX)
                    +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELX_INDEX);
                }
                else if(order_scheme_directional[0]==3)
                {
                	p.rdata(realData::xvel) = (alpha_pic_flip)*p.rdata(realData::xvel)
                	+(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELX_INDEX,lo,hi)
                	+(1-alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELX_INDEX,lo,hi);
                }

                if(order_scheme_directional[1]==1)
                {
                    p.rdata(realData::yvel) = (alpha_pic_flip)*p.rdata(realData::yvel)
                    +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELY_INDEX)
                    +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELY_INDEX);
                }
                else if(order_scheme_directional[1]==3)
                {
                	p.rdata(realData::yvel) = (alpha_pic_flip)*p.rdata(realData::yvel)
					+(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELY_INDEX,lo,hi)
                	+(1-alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELY_INDEX,lo,hi);
                }

                if(order_scheme_directional[2]==1)
                {
                    p.rdata(realData::zvel) = (alpha_pic_flip)*p.rdata(realData::zvel)
                    +(alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,DELTA_VELZ_INDEX)
                    +(1-alpha_pic_flip)*bilin_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELZ_INDEX);
                }
                if(order_scheme_directional[2]==3)
                {
                    p.rdata(realData::zvel) = (alpha_pic_flip)*p.rdata(realData::zvel)
                    +(alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,DELTA_VELZ_INDEX,lo,hi)
                    +(1-alpha_pic_flip)*cubic_interp(xp,iv[XDIR],iv[YDIR],iv[ZDIR],lmin,mmin,nmin,lmax,mmax,nmax,plo,dx,nodal_data_arr,VELZ_INDEX,lo,hi);
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

                int ind=0;
                for(int d1=0;d1<AMREX_SPACEDIM;d1++)
                {
                    for(int d2=d1;d2<AMREX_SPACEDIM;d2++)
                    {
                        p.rdata(realData::strainrate+ind)=0.5*(gradvp[d1][d2]+gradvp[d2][d1]);
                        /*if(ind==0)			//Use this block for gradient calculation checks
                        {
                        	p.rdata(realData::strainrate+ind)-=2.0*pi*cos(2.0*pi*p.pos(0));
                        }
                        if(ind==3)
                        {
                        	p.rdata(realData::strainrate+ind)-=2.0*pi*cos(2.0*pi*p.pos(1));
                        }
                        if(ind==5)
                        {
                        	p.rdata(realData::strainrate+ind)-=2.0*pi*cos(2.0*pi*p.pos(2));
                        }*/
                        p.rdata(realData::spinrate+ind)=0.5*(gradvp[d1][d2]-gradvp[d2][d1]); // this only calculates the upper half of spin tensor
                        ind++;
                    }
                }
            }
        });
    }
}
