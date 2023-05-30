#include <mpm_particle_container.H>
#include <interpolants.H>
#include <mpm_eb.H>
#include <mpm_kernels.H>
#include <constants.H>

amrex::Real MPMParticleContainer::Calculate_time_step(amrex::Real CFL,amrex::Real dtmax,amrex::Real dtmin)
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
        if(p.idata(intData::phase)==0)
        {
            amrex::Real velmag=std::sqrt(p.rdata(realData::xvel)*p.rdata(realData::xvel) +
										 p.rdata(realData::yvel)*p.rdata(realData::yvel) +
										 p.rdata(realData::zvel)*p.rdata(realData::zvel) );
			
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
            else if(p.idata(intData::constitutive_model)==2){
                if(velmag < TINYVAL){
                    Cs = TINYVAL;
                }
                else{
                    Cs = 0.0;
                }
            }

			Real tscale=amrex::min<amrex::Real>(dx[0],amrex::min<amrex::Real>(dx[1],dx[2]))/(Cs+velmag);
			return(tscale);
        }
        else
        {
        	Real tscale=std::numeric_limits<amrex::Real>::max();
        	return(tscale);
        }

    });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealMin(dt);
#endif

    dt	= CFL*dt;
    dt  = (dt>dtmax)?dtmax:((dt<dtmin)?dtmin:dt);


    return(dt);
}

void MPMParticleContainer::updateVolume(const amrex::Real& dt)
{
    BL_PROFILE("MPMParticleContainer::updateVolume");

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
            //yli add for gbhypo
            amrex::Real old_void_ratio = 0.0;
            amrex::Real new_void_ratio = 0.0;
            //end yli add
            if(p.idata(intData::phase)==0)
            {
                //yli add for gbhypo
                if(p.idata(intData::constitutive_model)==2){
                    old_void_ratio = p.rdata(realData::void_ratio);
                    new_void_ratio = old_void_ratio + (p.rdata(realData::strainrate+XX)
                                                    +p.rdata(realData::strainrate+YY)
                                                    +p.rdata(realData::strainrate+ZZ))
                                                    * (1 +old_void_ratio)* dt;
                    p.rdata(realData::void_ratio) = new_void_ratio;
                    p.rdata(realData::volume)	= (1+new_void_ratio)/(1+old_void_ratio)*p.rdata(realData::vol_init); // use order_of_scheme=1 to update vol_init
                    p.rdata(realData::density)	= p.rdata(realData::mass)/p.rdata(realData::volume);
                }
                //end yli add
                else{
                    p.rdata(realData::jacobian) = p.rdata(realData::deformation_gradient+0)*(p.rdata(realData::deformation_gradient+4)*p.rdata(realData::deformation_gradient+8)-p.rdata(realData::deformation_gradient+7)*p.rdata(realData::deformation_gradient+5))-
											  p.rdata(realData::deformation_gradient+1)*(p.rdata(realData::deformation_gradient+3)*p.rdata(realData::deformation_gradient+8)-p.rdata(realData::deformation_gradient+6)*p.rdata(realData::deformation_gradient+5))+
											  p.rdata(realData::deformation_gradient+2)*(p.rdata(realData::deformation_gradient+3)*p.rdata(realData::deformation_gradient+7)-p.rdata(realData::deformation_gradient+6)*p.rdata(realData::deformation_gradient+4));
                    p.rdata(realData::volume)	= p.rdata(realData::vol_init)*p.rdata(realData::jacobian);
                    p.rdata(realData::density)	= p.rdata(realData::mass)/p.rdata(realData::volume);
                }
            }
        });
    }
}

void MPMParticleContainer::moveParticles(const amrex::Real& dt,
        int bclo[AMREX_SPACEDIM],
		int bchi[AMREX_SPACEDIM],
		int lsetbc,
        amrex::Real wall_mu_lo[AMREX_SPACEDIM],
        amrex::Real wall_mu_hi[AMREX_SPACEDIM],
        amrex::Real wall_vel_lo[AMREX_SPACEDIM*AMREX_SPACEDIM],
        amrex::Real wall_vel_hi[AMREX_SPACEDIM*AMREX_SPACEDIM],
        amrex::Real lset_wall_mu)
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

    //Some GPU stuff again (Sreejith)
    GpuArray<int,AMREX_SPACEDIM> bc_lo_arr;
    GpuArray<int,AMREX_SPACEDIM> bc_hi_arr;

    for(int d=0;d<AMREX_SPACEDIM;d++)
    {
    	bc_lo_arr[d]=bclo[d];
    	bc_hi_arr[d]=bchi[d];
    }

    GpuArray<Real,AMREX_SPACEDIM*AMREX_SPACEDIM> wall_vel_lo_arr;
    GpuArray<Real,AMREX_SPACEDIM*AMREX_SPACEDIM> wall_vel_hi_arr;

    for(int d=0;d<AMREX_SPACEDIM*AMREX_SPACEDIM;d++)
    {
    	wall_vel_lo_arr[d]=wall_vel_lo[d];
    	wall_vel_hi_arr[d]=wall_vel_hi[d];
    }

    GpuArray<Real,AMREX_SPACEDIM*AMREX_SPACEDIM> wall_mu_lo_arr;
    GpuArray<Real,AMREX_SPACEDIM*AMREX_SPACEDIM> wall_mu_hi_arr;

    for(int d=0;d<AMREX_SPACEDIM*AMREX_SPACEDIM;d++)
    {
    	wall_mu_lo_arr[d]=wall_mu_lo[d];
    	wall_mu_hi_arr[d]=wall_mu_hi[d];
    }

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

            if(p.idata(intData::phase)==1)
            {
            	p.pos(XDIR) += p.rdata(realData::xvel_prime) * dt;
            	p.pos(YDIR) += p.rdata(realData::yvel_prime) * dt;
            	p.pos(ZDIR) += p.rdata(realData::zvel_prime) * dt;
            }

            if(p.idata(intData::phase)==0)
            {

            p.pos(XDIR) += p.rdata(realData::xvel_prime) * dt;
            p.pos(YDIR) += p.rdata(realData::yvel_prime) * dt;
            p.pos(ZDIR) += p.rdata(realData::zvel_prime) * dt;

            //for imposing boundary conditions           
            Real relvel_in[AMREX_SPACEDIM]  = {p.rdata(realData::xvel),p.rdata(realData::yvel),p.rdata(realData::zvel)};
            Real relvel_out[AMREX_SPACEDIM] = {p.rdata(realData::xvel),p.rdata(realData::yvel),p.rdata(realData::zvel)};

            if(using_levsets)
            {
                amrex::Real eps=0.00001;
                amrex::Real xp[AMREX_SPACEDIM]={p.pos(XDIR),p.pos(YDIR),p.pos(ZDIR)}; 
                amrex::Real dist=get_levelset_value(lsetarr,plo,dx,xp,lsref); 

               if(dist<TINYVAL)
               {
                    amrex::Real normaldir[AMREX_SPACEDIM]={1.0,0.0,0.0};
                    get_levelset_grad(lsetarr,plo,dx,xp,lsref,normaldir);
                    amrex::Real gradmag=std::sqrt(normaldir[XDIR]*normaldir[XDIR]
                                                 +normaldir[YDIR]*normaldir[YDIR]
                                                 +normaldir[ZDIR]*normaldir[ZDIR]);

                    normaldir[XDIR]=normaldir[XDIR]/(gradmag+TINYVAL);
                    normaldir[YDIR]=normaldir[YDIR]/(gradmag+TINYVAL);
                    normaldir[ZDIR]=normaldir[ZDIR]/(gradmag+TINYVAL);
                
                    int modify_pos=applybc(relvel_in,relvel_out,lset_wall_mu,
                        normaldir,lsetbc);
                    
                    if(modify_pos)
                    {
                        p.pos(XDIR) += 2.0*amrex::Math::abs(dist)*normaldir[XDIR];
                        p.pos(YDIR) += 2.0*amrex::Math::abs(dist)*normaldir[YDIR];
                        p.pos(ZDIR) += 2.0*amrex::Math::abs(dist)*normaldir[ZDIR];
                    }

                    p.rdata(realData::xvel)=relvel_out[XDIR];
                    p.rdata(realData::yvel)=relvel_out[YDIR];
                    p.rdata(realData::zvel)=relvel_out[ZDIR];
               }
            }
            
            relvel_in[XDIR]   = p.rdata(realData::xvel);
            relvel_in[YDIR]   = p.rdata(realData::yvel);
            relvel_in[ZDIR]   = p.rdata(realData::zvel);

            relvel_out[XDIR]  = p.rdata(realData::xvel);
            relvel_out[YDIR]  = p.rdata(realData::yvel);
            relvel_out[ZDIR]  = p.rdata(realData::zvel);

            amrex::Real wallvel[AMREX_SPACEDIM]={0.0,0.0,0.0};

            if (p.pos(XDIR) < plo[XDIR])
            {
                int dir=XDIR;
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    wallvel[d]=wall_vel_lo_arr[dir*AMREX_SPACEDIM+d];
                    relvel_in[d] -= wallvel[d];
                }
                
                Real normaldir[AMREX_SPACEDIM]={1.0,0.0,0.0};
                //int modify_pos=applybc(relvel_in,relvel_out,wall_mu_lo[XDIR],normaldir,bclo[XDIR]);
                int modify_pos=applybc(relvel_in,relvel_out,wall_mu_lo_arr[XDIR],normaldir,bc_lo_arr[XDIR]);
                if(modify_pos)
                {
                    p.pos(XDIR) = two*plo[XDIR] - p.pos(XDIR);
                }
            }
            else if (p.pos(XDIR) > phi[XDIR])
            {
                int dir=XDIR;
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    wallvel[d]=wall_vel_hi_arr[dir*AMREX_SPACEDIM+d];
                    relvel_in[d] -= wallvel[d];
                }

                Real normaldir[AMREX_SPACEDIM]={-1.0,0.0,0.0};
                //int modify_pos=applybc(relvel_in,relvel_out,wall_mu_hi[XDIR],normaldir,bchi[XDIR]);
                int modify_pos=applybc(relvel_in,relvel_out,wall_mu_hi_arr[XDIR],normaldir,bc_hi_arr[XDIR]);
                if(modify_pos)
                {
                    p.pos(XDIR) = two*phi[XDIR] - p.pos(XDIR);
                }
            }
            if (p.pos(YDIR) < plo[YDIR])
            {
                int dir=YDIR;
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    wallvel[d]=wall_vel_lo_arr[dir*AMREX_SPACEDIM+d];
                    relvel_in[d] -= wallvel[d];
                }
                
                Real normaldir[AMREX_SPACEDIM]={0.0,1.0,0.0};
                //int modify_pos=applybc(relvel_in,relvel_out,wall_mu_lo[YDIR],normaldir,bclo[YDIR]);
                int modify_pos=applybc(relvel_in,relvel_out,wall_mu_lo_arr[YDIR],normaldir,bc_lo_arr[YDIR]);
                if(modify_pos)
                {
                    p.pos(YDIR) = two*plo[YDIR] - p.pos(YDIR);
                }
            }
            else if (p.pos(YDIR) > phi[YDIR])
            {
                int dir=YDIR;
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    wallvel[d]=wall_vel_hi_arr[dir*AMREX_SPACEDIM+d];
                    relvel_in[d] -= wallvel[d];
                }
            
                Real normaldir[AMREX_SPACEDIM]={0.0,-1.0,0.0};
                //int modify_pos=applybc(relvel_in,relvel_out,wall_mu_hi[YDIR],normaldir,bchi[YDIR]);
                int modify_pos=applybc(relvel_in,relvel_out,wall_mu_hi[YDIR],normaldir,bc_hi_arr[YDIR]);
                if(modify_pos)
                {
                    p.pos(YDIR) = two*phi[YDIR] - p.pos(YDIR);
                }
            }
            if (p.pos(ZDIR) < plo[ZDIR])
            {
                int dir=ZDIR;
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    wallvel[d]=wall_vel_lo_arr[dir*AMREX_SPACEDIM+d];
                    relvel_in[d] -= wallvel[d];
                }
            
                Real normaldir[AMREX_SPACEDIM]={0.0,0.0,1.0};
                //int modify_pos=applybc(relvel_in,relvel_out,wall_mu_lo[ZDIR],normaldir,bclo[ZDIR]);
                int modify_pos=applybc(relvel_in,relvel_out,wall_mu_lo_arr[ZDIR],normaldir,bc_lo_arr[ZDIR]);
                if(modify_pos)
                {
                    p.pos(ZDIR) = two*plo[ZDIR] - p.pos(ZDIR);
                }
            }
            else if (p.pos(ZDIR) > phi[ZDIR])
            {
                int dir=ZDIR;
                for(int d=0;d<AMREX_SPACEDIM;d++)
                {
                    wallvel[d]=wall_vel_hi_arr[dir*AMREX_SPACEDIM+d];
                    relvel_in[d] -= wallvel[d];
                }
            
                Real normaldir[AMREX_SPACEDIM]={0.0,0.0,-1.0};
                //int modify_pos=applybc(relvel_in,relvel_out,wall_mu_hi[ZDIR],normaldir,bchi[ZDIR]);
                int modify_pos=applybc(relvel_in,relvel_out,wall_mu_hi_arr[ZDIR],normaldir,bc_hi_arr[ZDIR]);
                if(modify_pos)
                {
                    p.pos(ZDIR) = two*phi[ZDIR] - p.pos(ZDIR);
                }
            }
            else //nothing to do
            {}
            p.rdata(realData::xvel)=relvel_out[XDIR]+wallvel[XDIR];
            p.rdata(realData::yvel)=relvel_out[YDIR]+wallvel[YDIR];
            p.rdata(realData::zvel)=relvel_out[ZDIR]+wallvel[ZDIR];
            }
        });
    }
}

amrex::Real MPMParticleContainer::GetPosSpring()
{
	const int lev = 0;
	const Geometry& geom = Geom(lev);
	const auto plo = Geom(lev).ProbLoArray();
	const auto phi = Geom(lev).ProbHiArray();
	const auto dx = Geom(lev).CellSizeArray();
	auto& plev  = GetParticles(lev);
	amrex::Real ymin = 0.0;

	using PType = typename MPMParticleContainer::SuperParticleType;
	ymin = amrex::ReduceMax(*this, [=]
			AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
	        {
	        	Real yscale;
	        	yscale = p.pos(YDIR);
	        	return(yscale);
	         });
	return(ymin);

}

amrex::Real MPMParticleContainer::GetPosPiston()
{
	const int lev = 0;
	const Geometry& geom = Geom(lev);
	const auto plo = Geom(lev).ProbLoArray();
	const auto phi = Geom(lev).ProbHiArray();
	const auto dx = Geom(lev).CellSizeArray();
	auto& plev  = GetParticles(lev);
	amrex::Real ymin = std::numeric_limits<amrex::Real>::max();

	using PType = typename MPMParticleContainer::SuperParticleType;
	ymin = amrex::ReduceMin(*this, [=]
			AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
	        {
	        	Real yscale;
	        	if(p.idata(intData::phase)==1 and p.idata(intData::rigid_body_id)==0)
	            {
	        		yscale = p.pos(YDIR);
	            }
	        	else
	        	{
	        		yscale = std::numeric_limits<amrex::Real>::max();
	        	}
	        	return(yscale);
	         });
	return(ymin);

}


void MPMParticleContainer::UpdateRigidParticleVelocities(int rigid_body_id,Array <amrex::Real,AMREX_SPACEDIM> velocity)
{
    BL_PROFILE("MPMParticleContainer::GetVelPiston");

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
                if(p.idata(intData::phase)==1 and p.idata(intData::rigid_body_id)==rigid_body_id)
                {
                	p.rdata(realData::xvel_prime) =velocity[0];
                	p.rdata(realData::yvel_prime) =velocity[1];
                	p.rdata(realData::zvel_prime) =velocity[2];
                }
            });
        }

}


