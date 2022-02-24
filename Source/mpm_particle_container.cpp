#include <mpm_particle_container.H>
#include <interpolants.H>
#include <constitutive_models.H>

using namespace amrex;

void MPMParticleContainer::apply_constitutive_model(const amrex::Real& dt,
                                                    amrex::Real E,amrex::Real v,
                                                    amrex::Real applied_strainrate=0.0)
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

            linear_elastic(strain,strainrate,stress,E,v);

            for(int d=0;d<NCOMP_TENSOR;d++)
            {
                p.rdata(realData::stress+d)=stress[d];
            }
        });
    }
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

void MPMParticleContainer::deposit_onto_grid(MultiFab& nodaldata,
                                             Array<Real,AMREX_SPACEDIM> gravity,
                                             int update_massvel,int update_forces) 
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();
    Real grav[]={AMREX_D_DECL(gravity[0],gravity[1],gravity[2])};

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
            ParticleType& p = pstruct[i];

            amrex::Real xp[AMREX_SPACEDIM];

            xp[XDIR]=p.pos(XDIR);
            xp[YDIR]=p.pos(YDIR);
            xp[ZDIR]=p.pos(ZDIR);

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

                            amrex::Real basisvalue=basisval(l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx);

                            if(update_massvel)
                            {

                                amrex::Real mass_contrib=p.rdata(realData::mass)*basisvalue;

                                amrex::Real p_contrib[AMREX_SPACEDIM] = 
                                {p.rdata(realData::mass)*p.rdata(realData::xvel)*basisvalue,
                                    p.rdata(realData::mass)*p.rdata(realData::yvel)*basisvalue,
                                    p.rdata(realData::mass)*p.rdata(realData::zvel)*basisvalue};

                                amrex::Gpu::Atomic::AddNoRet(
                                    &nodal_data_arr(ivlocal,MASS_INDEX), mass_contrib);

                                for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                                {
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &nodal_data_arr(ivlocal,VELX_INDEX+dim), 
                                        p_contrib[dim]);
                                }
                            }

                            if(update_forces)
                            {
                                amrex::Real basisval_grad[AMREX_SPACEDIM];
                                amrex::Real stress_tens[AMREX_SPACEDIM*AMREX_SPACEDIM];

                                get_tensor(p,realData::stress,stress_tens);

                                for(int d=0;d<AMREX_SPACEDIM;d++)
                                {
                                    basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx);
                                }

                                amrex::Real bforce_contrib[AMREX_SPACEDIM]=
                                {p.rdata(realData::mass)*grav[XDIR]*basisvalue,
                                    p.rdata(realData::mass)*grav[YDIR]*basisvalue,
                                    p.rdata(realData::mass)*grav[ZDIR]*basisvalue};

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

        amrex::ParallelFor(
            nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if(update_massvel)
                {
                    if(nodal_data_arr(i,j,k,MASS_INDEX) > 0.0)
                    {
                        for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                        {
                            nodal_data_arr(i,j,k,VELX_INDEX+dim)/=nodal_data_arr(i,j,k,MASS_INDEX);
                        }
                    }
                }
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

            if ( p.pos(0) < plo[0])
            {
                p.pos(0) = two*plo[0] - p.pos(0);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (p.pos(0) > phi[0])
            {
                p.pos(0) = two*phi[0] - p.pos(0);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (p.pos(1) < plo[1])
            {
                p.pos(1) = two*plo[1] - p.pos(1);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (p.pos(1) > phi[1])
            {
                p.pos(1) = two*phi[1] - p.pos(1);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (p.pos(2) < plo[2])
            {
                p.pos(2) = two*plo[2] - p.pos(2);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }
            if (p.pos(2) > phi[2])
            {
                p.pos(2) = two*phi[2] - p.pos(2);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }

        });
    }
}

void MPMParticleContainer::interpolate_from_grid(MultiFab& nodaldata,int update_vel,int update_strainrate)
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

            if(update_vel)
            {
                p.rdata(realData::xvel) = bilin_interp(xp,iv[XDIR],
                             iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELX_INDEX);

                p.rdata(realData::yvel) = bilin_interp(xp,iv[XDIR],
                             iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELY_INDEX);

                p.rdata(realData::zvel) = bilin_interp(xp,iv[XDIR],
                             iv[YDIR],iv[ZDIR],plo,dx,nodal_data_arr,VELZ_INDEX);
            }

            if(update_strainrate)
            {
                for(int n=0;n<2;n++)
                {
                    for(int m=0;m<2;m++)
                    {
                        for(int l=0;l<2;l++)
                        {
                            amrex::Real basisval_grad[AMREX_SPACEDIM];
                            for(int d=0;d<AMREX_SPACEDIM;d++)
                            {
                                basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx);
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

    int_data_names.push_back("phase");

    writeflags_int[intData::phase]=1;

    writeflags_real[realData::radius]=1;
    writeflags_real[realData::xvel]=1;
    writeflags_real[realData::yvel]=1;
    writeflags_real[realData::zvel]=1;
    writeflags_real[realData::mass]=1;

    WritePlotFile(pltfile, "particles",writeflags_real, 
                  writeflags_int, real_data_names, int_data_names);
}
