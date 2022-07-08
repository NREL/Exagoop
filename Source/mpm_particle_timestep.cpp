#include <mpm_particle_container.H>
#include <interpolants.H>
#include <mpm_eb.H>

amrex::Real MPMParticleContainer::Calculate_time_step()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dx = geom.CellSizeArray();
    amrex::Real dt = std::numeric_limits<amrex::Real>::max();
    
    using PType = typename MPMParticleContainer::SuperParticleType;
    dt = amrex::ReduceMin(*this, [=] 
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real 
    {
        amrex::Real Cs;
        if(p.idata(intData::constitutive_model)==1)
        {
            Cs = sqrt(p.rdata(realData::Bulk_modulus)/p.rdata(realData::density));
        }
        else if(p.idata(intData::constitutive_model)==0)
        {
            Real lambda=p.rdata(realData::E)*p.rdata(realData::nu)/
            ((1+p.rdata(realData::nu))*(1-2.0*p.rdata(realData::nu)));
            Real mu=p.rdata(realData::E)/(2.0*(1+p.rdata(realData::nu)));
            Cs = sqrt((lambda+2.0*mu)/p.rdata(realData::density));
        }
        amrex::Real velmag=std::sqrt(p.rdata(realData::xvel)*p.rdata(realData::xvel) +
                                     p.rdata(realData::yvel)*p.rdata(realData::yvel) +
                                     p.rdata(realData::zvel)*p.rdata(realData::zvel) );

        Real tscale=amrex::min<amrex::Real>(dx[0],amrex::min<amrex::Real>(dx[1],dx[2]))/(Cs+velmag);
        return(tscale);
    });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealMin(dt);
#endif

    if(dt<1e-10)
    {
        amrex::Print()<<"\nWarning: Time step is getting too low (dt = "<<dt<<" )";
    }
    return(dt);
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

            p.rdata(realData::jacobian) += (p.rdata(realData::strainrate+XX)
            +p.rdata(realData::strainrate+YY)+p.rdata(realData::strainrate+ZZ)) * dt * p.rdata(realData::jacobian);
            p.rdata(realData::volume)	= p.rdata(realData::vol_init)*p.rdata(realData::jacobian);
            p.rdata(realData::density)	= p.rdata(realData::mass)/p.rdata(realData::volume);
            // Yudong: add void ratio evolution for hypoplastic model
            p.rdata(realData::void_ratio) += (p.rdata(realData::strainrate+XX)
            +p.rdata(realData::strainrate+YY)+p.rdata(realData::strainrate+ZZ))  * (1+p.rdata(realData::void_ratio))* dt;

        });
    }
}

void MPMParticleContainer::moveParticles(const amrex::Real& dt,
        int bclo[AMREX_SPACEDIM],int bchi[AMREX_SPACEDIM])
{
    BL_PROFILE("MPMParticleContainer::moveParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
    auto& plev  = GetParticles(lev);

    bool using_levsets=mpm_ebtools::using_levelset_geometry;
    int lsref=mpm_ebtools::ls_refinement;

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

        amrex::Array4<amrex::Real> lsetarr;
        if(using_levsets)
        {
            lsetarr=mpm_ebtools::lsphi->array(mfi);
        }

        // now we move the particles
        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];

            p.pos(XDIR) += p.rdata(realData::xvel) * dt;
            p.pos(YDIR) += p.rdata(realData::yvel) * dt;
            p.pos(ZDIR) += p.rdata(realData::zvel) * dt;

            if(using_levsets)
            {
                amrex::Real xp[AMREX_SPACEDIM]={p.pos(XDIR),p.pos(YDIR),p.pos(ZDIR)}; 
                amrex::Real dist=get_levelset_value(lsetarr,plo,dx,xp,lsref); 

               if(dist<TINYVAL)
               {
                    amrex::Real norm[AMREX_SPACEDIM]={1.0,0.0,0.0};
                    get_levelset_grad(lsetarr,plo,dx,xp,lsref,norm);

                    amrex::Real gradmag=std::sqrt(norm[XDIR]*norm[XDIR]+norm[YDIR]*norm[YDIR]+norm[ZDIR]*norm[ZDIR]);
                    amrex::Real velmag=std::sqrt(p.rdata(realData::xvel)*p.rdata(realData::xvel) +
                                     p.rdata(realData::yvel)*p.rdata(realData::yvel) +
                                     p.rdata(realData::zvel)*p.rdata(realData::zvel));

                    if(gradmag > TINYVAL)
                    {
                        norm[XDIR]=norm[XDIR]/gradmag;
                        norm[YDIR]=norm[YDIR]/gradmag;
                        norm[ZDIR]=norm[ZDIR]/gradmag;
                    }
                    p.pos(XDIR)+=2.0*amrex::Math::abs(dist)*norm[XDIR];
                    p.pos(YDIR)+=2.0*amrex::Math::abs(dist)*norm[YDIR];
                    p.pos(ZDIR)+=2.0*amrex::Math::abs(dist)*norm[ZDIR];

                    p.rdata(realData::xvel)=velmag*norm[XDIR];
                    p.rdata(realData::yvel)=velmag*norm[YDIR];
                    p.rdata(realData::zvel)=velmag*norm[ZDIR];
               }
            }

            if (!periodic[XDIR] && (p.pos(XDIR) < plo[XDIR]) && bclo[XDIR]!=BC_OUTFLOW)
            {
                p.pos(XDIR) = two*plo[XDIR] - p.pos(XDIR);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (!periodic[XDIR] && (p.pos(XDIR) > phi[XDIR]) && bchi[XDIR]!=BC_OUTFLOW)
            {
                p.pos(XDIR) = two*phi[XDIR] - p.pos(XDIR);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (!periodic[YDIR] && (p.pos(YDIR) < plo[YDIR]) && bclo[YDIR]!=BC_OUTFLOW)
            {
                p.pos(YDIR) = two*plo[YDIR] - p.pos(YDIR);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (!periodic[YDIR] && (p.pos(YDIR) > phi[YDIR]) && bchi[YDIR]!=BC_OUTFLOW)
            {
                p.pos(YDIR) = two*phi[YDIR] - p.pos(YDIR);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (!periodic[ZDIR] && (p.pos(ZDIR) < plo[ZDIR]) && bclo[ZDIR]!=BC_OUTFLOW)
            {
                p.pos(ZDIR) = two*plo[ZDIR] - p.pos(ZDIR);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }
            if (!periodic[ZDIR] && (p.pos(ZDIR) > phi[ZDIR]) && bchi[ZDIR]!=BC_OUTFLOW)
            {
                p.pos(ZDIR) = two*phi[ZDIR] - p.pos(ZDIR);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }

        });
    }
}

void MPMParticleContainer::move_particles_from_nodevel(MultiFab& nodaldata,
        const amrex::Real& dt,
        int bclo[AMREX_SPACEDIM],int bchi[AMREX_SPACEDIM],
        int order_scheme)
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
            
            if (!periodic[XDIR] && (p.pos(XDIR) < plo[XDIR]) && bclo[XDIR]!=BC_OUTFLOW)
            {
                p.pos(XDIR) = two*plo[XDIR] - p.pos(XDIR);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (!periodic[XDIR] && (p.pos(XDIR) > phi[XDIR]) && bchi[XDIR]!=BC_OUTFLOW)
            {
                p.pos(XDIR) = two*phi[XDIR] - p.pos(XDIR);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (!periodic[YDIR] && (p.pos(YDIR) < plo[YDIR]) && bclo[YDIR]!=BC_OUTFLOW)
            {
                p.pos(YDIR) = two*plo[YDIR] - p.pos(YDIR);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (!periodic[YDIR] && (p.pos(YDIR) > phi[YDIR]) && bchi[YDIR]!=BC_OUTFLOW)
            {
                p.pos(YDIR) = two*phi[YDIR] - p.pos(YDIR);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (!periodic[ZDIR] && (p.pos(ZDIR) < plo[ZDIR]) && bclo[ZDIR]!=BC_OUTFLOW)
            {
                p.pos(ZDIR) = two*plo[ZDIR] - p.pos(ZDIR);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }
            if (!periodic[ZDIR] && (p.pos(ZDIR) > phi[ZDIR]) && bchi[ZDIR]!=BC_OUTFLOW)
            {
                p.pos(ZDIR) = two*phi[ZDIR] - p.pos(ZDIR);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }

        });
    }
}
