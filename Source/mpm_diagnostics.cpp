#include <mpm_particle_container.H>
#include <interpolants.H>


void MPMParticleContainer::CalculateSurfaceLoads(Real &Load,int rigidbodyid)
{
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  auto& plev  = GetParticles(lev);
  const auto dxi = geom.InvCellSizeArray();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto domain = geom.Domain();

  Load=0.0;

  using PType = typename MPMParticleContainer::SuperParticleType;
  Load = amrex::ReduceSum(*this, [=]
                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                  {
    if(p.idata(intData::phase)==1 and  p.idata(intData::rigid_body_id)==rigidbodyid)
      {
        return(p.rdata(realData::stress+YY));
      }
    else
      {
        return(0.0);
      }
                                  });

#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealSum(Load);
#endif

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

void MPMParticleContainer::CalculateStressDiagnostics(amrex::GpuArray<amrex::Real, 6> &min_stress, amrex::GpuArray<amrex::Real, 6> &max_stress, amrex::GpuArray<amrex::Real, 6> &avg_stress)
{
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  auto& plev  = GetParticles(lev);
  const auto dxi = geom.InvCellSizeArray();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto domain = geom.Domain();

  using PType = typename MPMParticleContainer::SuperParticleType;
  for(int comp=0;comp<NCOMP_TENSOR;comp++)
    {
      min_stress[comp] = amrex::ReduceMin(*this, [=]
                                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                                  {
        return(p.rdata(realData::stress+comp));
                                                  });

    }
  for(int comp=0;comp<NCOMP_TENSOR;comp++)
    {
      max_stress[comp] = amrex::ReduceMax(*this, [=]
                                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                                  {
        return(p.rdata(realData::stress+comp));
                                                  });
    }
  for(int comp=0;comp<NCOMP_TENSOR;comp++)
    {
      avg_stress[comp] = amrex::ReduceSum(*this, [=]
                                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                                  {
        return(p.rdata(realData::stress+comp));
                                                  });
    }

#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealMin(min_stress[0]);
  ParallelDescriptor::ReduceRealMin(min_stress[1]);
  ParallelDescriptor::ReduceRealMin(min_stress[2]);
  ParallelDescriptor::ReduceRealMin(min_stress[3]);
  ParallelDescriptor::ReduceRealMin(min_stress[4]);
  ParallelDescriptor::ReduceRealMin(min_stress[5]);

  ParallelDescriptor::ReduceRealMax(max_stress[0]);
  ParallelDescriptor::ReduceRealMax(max_stress[1]);
  ParallelDescriptor::ReduceRealMax(max_stress[2]);
  ParallelDescriptor::ReduceRealMax(max_stress[3]);
  ParallelDescriptor::ReduceRealMax(max_stress[4]);
  ParallelDescriptor::ReduceRealMax(max_stress[5]);

  ParallelDescriptor::ReduceRealSum(avg_stress[0]);
  ParallelDescriptor::ReduceRealSum(avg_stress[1]);
  ParallelDescriptor::ReduceRealSum(avg_stress[2]);
  ParallelDescriptor::ReduceRealSum(avg_stress[3]);
  ParallelDescriptor::ReduceRealSum(avg_stress[4]);
  ParallelDescriptor::ReduceRealSum(avg_stress[5]);

#endif
}


void MPMParticleContainer::CalculatePosDiagnostics(amrex::GpuArray<amrex::Real, 3> &min_pos, amrex::GpuArray<amrex::Real, 3> &max_pos)
{
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  auto& plev  = GetParticles(lev);
  const auto dxi = geom.InvCellSizeArray();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto domain = geom.Domain();

  using PType = typename MPMParticleContainer::SuperParticleType;
  for(int comp=0;comp<AMREX_SPACEDIM;comp++)
    {
      min_pos[comp] = amrex::ReduceMin(*this, [=]
                                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                                  {
        return(p.pos(comp));
                                                  });

    }
  for(int comp=0;comp<AMREX_SPACEDIM;comp++)
    {
      max_pos[comp] = amrex::ReduceMax(*this, [=]
                                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                                  {
        return(p.pos(comp));
                                                  });
    }

#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealMin(min_pos[0]);
  ParallelDescriptor::ReduceRealMin(min_pos[1]);
  ParallelDescriptor::ReduceRealMin(min_pos[2]);

  ParallelDescriptor::ReduceRealMax(max_pos[0]);
  ParallelDescriptor::ReduceRealMax(max_pos[1]);
  ParallelDescriptor::ReduceRealMax(max_pos[2]);

#endif
}

void MPMParticleContainer::CalculateVelocityDiagnostics(Real &min_vel,Real &max_vel, Real &avg_vel)
{
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  auto& plev  = GetParticles(lev);
  const auto dxi = geom.InvCellSizeArray();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto domain = geom.Domain();

  min_vel=0.0;
  max_vel=0.0;
  avg_vel=0.0;

  using PType = typename MPMParticleContainer::SuperParticleType;
  min_vel = amrex::ReduceMin(*this, [=]
                                     AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                     {
    return(sqrt
        (p.rdata(realData::xvel)*p.rdata(realData::xvel)+
         p.rdata(realData::yvel)*p.rdata(realData::yvel)+
         p.rdata(realData::zvel)*p.rdata(realData::zvel)) );
                                     });

  max_vel = amrex::ReduceMax(*this, [=]
                                     AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                     {
    return(sqrt
        (p.rdata(realData::xvel)*p.rdata(realData::xvel)+
         p.rdata(realData::yvel)*p.rdata(realData::yvel)+
         p.rdata(realData::zvel)*p.rdata(realData::zvel)) );
                                     });

  avg_vel = amrex::ReduceSum(*this, [=]                                       //This is actually the sum. Divide it by the total number of material points in the main function
                                     AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                     {
    return(sqrt
        (p.rdata(realData::xvel)*p.rdata(realData::xvel)+
         p.rdata(realData::yvel)*p.rdata(realData::yvel)+
         p.rdata(realData::zvel)*p.rdata(realData::zvel)) );
                                     });

#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealMin(min_vel);
  ParallelDescriptor::ReduceRealMax(max_vel);
  ParallelDescriptor::ReduceRealSum(avg_vel);
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

void MPMParticleContainer::ParticleErrorAxialBar(MultiFab& nodaldata,
                                                 int n, Real L, Real E, Real rho, Real time, Real V0)
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
    const double Pi = 3.141592654;

    Real beta_n=(2.0*n-1)/2.0*Pi/L;
    Real omega_n=beta_n*sqrt(E/rho);



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
                Real Vex=0.0;
                Vex=V0*sin(Pi*p.pos(0)/(2.0*L))*cos(omega_n*time);
                p.rdata(realData::error) = fabs(p.rdata(realData::xvel)-Vex);
            }
        });
    }
}

void MPMParticleContainer::ParticleErrorTranslationofNonInteractingMaterialPoints(MultiFab& nodaldata,
                                                 Real V0, Real g,Real time)
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
    const double Pi = 3.141592654;





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
                Real position=0.0;
                position=p.rdata(realData::pos_0_y)+V0*time+0.5*g*time*time;
                p.rdata(realData::error) = fabs(p.pos(1)-position);
            }
        });
    }
}


void MPMParticleContainer::CalculatePositionL2Error_TranslationofNonInteractingMaterialPoints(Real &Error, Real V0, Real G0, Real W0, Real time, int nump)
{
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  auto& plev  = GetParticles(lev);
  const auto dxi = geom.InvCellSizeArray();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto domain = geom.Domain();

  Real Pi = atan(1.0)*4.0;



  using PType = typename MPMParticleContainer::SuperParticleType;
  Error = amrex::ReduceSum(*this, [=]
                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                  {
    Real position=0.0;
    if(W0==0)
      {
        position= p.rdata(realData::pos_0_y)+V0*time+0.5*G0*time*time;               //s=x0+u*t+1/2*g*t^2. Only y-displacement is computed.
      }
    else
      {
        position= p.rdata(realData::pos_0_y)+G0/(W0*W0)+V0*time-G0/(W0*W0)*cos(W0*time);               //s=x0+u*t+1/2*g*t^2. Only y-displacement is computed.
      }

    return((p.pos(1)-position)*(p.pos(1)-position));
                                  });



#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealSum(Error);
#endif

  Error=sqrt(Error)/nump;


}

void MPMParticleContainer::CalculateVelocityL2Error_Axialbar(Real &Error, int n, Real L, Real E, Real rho, Real time, Real V0, int nump)
{
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  auto& plev  = GetParticles(lev);
  const auto dxi = geom.InvCellSizeArray();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto domain = geom.Domain();

  Real Pi = atan(1.0)*4.0;
  Real beta_n=(2.0*n-1)/2.0*Pi/L;
  Real omega_n=beta_n*sqrt(E/rho);


  using PType = typename MPMParticleContainer::SuperParticleType;
  Error = amrex::ReduceSum(*this, [=]
                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                  {
    Real Vex=0.0;
    Vex=V0*sin(Pi*p.pos(0)/(2.0*L))*cos(omega_n*time);
    //p.rdata(realData::error) = p.rdata(realData::xvel)-Vex;
    return((p.rdata(realData::xvel)-Vex)*(p.rdata(realData::xvel)-Vex));
                                  });



#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealSum(Error);
#endif

  Error=sqrt(Error)/nump;
}

void MPMParticleContainer::CalculateVelocityL2Error_TransverseVibrationString(Real &Error, Real L, Real time, Real V0, Real c, int nump)
{
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  auto& plev  = GetParticles(lev);
  const auto dxi = geom.InvCellSizeArray();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto domain = geom.Domain();

  Real Pi = atan(1.0)*4.0;

  using PType = typename MPMParticleContainer::SuperParticleType;
  Error = amrex::ReduceSum(*this, [=]
                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                  {
    Real Yex=V0*L/(c*Pi)*sin(c*Pi*time/L)*sin(Pi*p.pos(0)/L);
    return((p.pos(1)-Yex)*(p.pos(1)-Yex));
                                  });



#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealSum(Error);
#endif

  Error=sqrt(Error)/nump;


}

void MPMParticleContainer::CalculatePositionL2Error_PositionofSingleMaterialPoint(Real &Error, Real x0, Real u0, Real accln, Real time, int nump)
{
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  auto& plev  = GetParticles(lev);
  const auto dxi = geom.InvCellSizeArray();
  const auto dx = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();
  const auto domain = geom.Domain();



  using PType = typename MPMParticleContainer::SuperParticleType;
  Error = amrex::ReduceSum(*this, [=]
                                  AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                  {
    Real Xex=x0+ u0*time+0.5*accln*time*time;
    return((p.pos(0)-Xex)*(p.pos(0)-Xex));
                                  });

#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealSum(Error);
#endif

  Error=sqrt(Error)/nump;
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
        //PrintToFile("CantileverDeflection.out")<<xp[XDIR]<<"\t"<<xp[YDIR]<<"\n";
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

  const int* loarr = domain.loVect ();
  const int* hiarr = domain.hiVect ();

  int lo[]={loarr[0],loarr[1],loarr[2]};
  int hi[]={hiarr[0],hiarr[1],hiarr[2]};

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
            //PrintToFile(outputfile)<<x<<"\t"<<v<<"\n";
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
        //PrintToFile(outputfile).SetPrecision(17)<<xp[XDIR]<<"\t"<<y_exact<<"\n";
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
  amrex::Print()<<"\n Exact = "<<tvb_v0<<" "<<tvb_L<<" "<<tvb_E<<" "<<tvb_rho<<" "<<pi<<" "<<w0;
  std::string outputfile;
  int num_of_digits_in_filenames = 5;
  outputfile = amrex::Concatenate("TVB_Deflection_", output_it, num_of_digits_in_filenames);
  //amrex::Print()<<"\nE = "<<tvb_E<<" "<<tvb_rho<<" "<<c<<" "<<w0;

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
#ifndef AMREX_USE_GPU
        PrintToFile(outputfile)<<xp[XDIR]<<"\t"<<xp[YDIR]<<"\t"<<y_exact<<"\n";
#endif
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

#ifndef AMREX_USE_GPU
  PrintToFile("SpringConst.out")<<Calculated_Spring_Const<<"\n";
#endif

  //Calculate and return the restoring force
  return(Restoring_force);

}


