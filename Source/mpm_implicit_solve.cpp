#include <mpm_particle_container.H>
#include <mpm_implicit.H>
#include <cvode/cvode.h>
#include <interpolants.H>
#include <constitutive_models.H>

void MPMParticleContainer::implicitUpdate(int implicit_solve, 
                                          Vector<MultiFab>& v_vect_old, 
                                          Vector<MultiFab>& v_vect_new, 
                                          MultiFab& nodaldata,
                                          Array<Real,AMREX_SPACEDIM> gravity,
                                          int external_loads_present,
                                          Array<Real,AMREX_SPACEDIM> force_slab_lo,
                                          Array<Real,AMREX_SPACEDIM> force_slab_hi,
                                          Array<Real,AMREX_SPACEDIM> extforce,
                                          Real mass_tolerance,
                                          GpuArray<int, AMREX_SPACEDIM> order_scheme_directional,
                                          GpuArray<int, AMREX_SPACEDIM> periodic,
                                          amrex::Real applied_strainrate,
                                          Real time, Real dt){

    if(implicit_solve == 1){
        implicitSolveM1(v_vect_old, v_vect_new, 
                        nodaldata,  
                        gravity,
                        external_loads_present,
                        force_slab_lo,
                        force_slab_hi,
                        extforce,
                        mass_tolerance, 
                        order_scheme_directional,
                        periodic, applied_strainrate,
                        time, dt);    
    } else if(implicit_solve == 2){
    
    } else {
        Error("ERROR: Invalid implicit method selected!\n");
    }
}

void MPMParticleContainer::implicitSolveM1(Vector<MultiFab>& v_vect_old, 
                                           Vector<MultiFab>& v_vect_new,
                                           MultiFab& nodaldata, 
                                           Array<Real,AMREX_SPACEDIM> gravity,
                                           int external_loads_present,
                                           Array<Real,AMREX_SPACEDIM> force_slab_lo,
                                           Array<Real,AMREX_SPACEDIM> force_slab_hi,
                                           Array<Real,AMREX_SPACEDIM> extforce,
                                           Real mass_tolerance, 
                                           GpuArray<int, AMREX_SPACEDIM> order_scheme_directional,
                                           GpuArray<int, AMREX_SPACEDIM> periodic,
                                           amrex::Real applied_strainrate,
                                           Real time, Real dt){

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

    // Define RHS function for the AMReX integrator
    auto rhs_function = [&](Vector<MultiFab>& S_rhs,
                      const Vector<MultiFab>& S_data, const Real /* time */) {

        auto& v_data = S_data[0];
        auto& v_rhs  = S_rhs[0];

        /* ------------------------------------------------------------*
         *                  Calculate internal forces                  *
         * ------------------------------------------------------------*/

        // First zero out any nodal variables used in RHS evaluation
        // TODO: Do we need to zero anything else?
        for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
        {
            const Box& nodalbox=mfi.validbox();
            Array4<Real> nodal_data_arr=nodaldata.array(mfi);

            amrex::ParallelFor(nodalbox,[=]
            AMREX_GPU_DEVICE (int i,int j,int k) noexcept
            {
                nodal_data_arr(i,j,k,FRCX_INDEX)=ZERO;
                nodal_data_arr(i,j,k,FRCY_INDEX)=ZERO;
                nodal_data_arr(i,j,k,FRCZ_INDEX)=ZERO;
            });
        }

        // Loop over real particles to update volume and stress
        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){
            const amrex::Box& box = mfi.tilebox();
            Box nodalbox = convert(box, {1, 1, 1});

            int gid = mfi.index();
            int tid = mfi.LocalTileIndex();
            auto index = std::make_pair(gid, tid);

            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            int np = aos.numRealParticles();
            // int ng =aos.numNeighborParticles();
            // int nt = np+ng;

            Array4<Real> nodal_data_arr=nodaldata.array(mfi);
            const Array4<const Real>& v_array = v_data.array(mfi);
            ParticleType* pstruct = aos().dataPtr();

            // Begin be looping over all grid points to update the 
            // particle volume and stress tensors
            // TODO: Do we loop over only real or real + neighbor particles here?
            //        Don't want to double count updates but also need to ensure necessary data is available..
            amrex::ParallelFor(np,[=]
            AMREX_GPU_DEVICE (int i) noexcept{ 
                int lmin,lmax,nmin,nmax,mmin,mmax;
                ParticleType& p = pstruct[i];
                if(p.idata(intData::phase)==0)    //Compute only for standard particles and not rigid particles with phase=1
                {
                    // Reset values used to calculate volume and stresses to t=n
                    for(int d=0; d<NCOMP_FULLTENSOR; d++) p.rdata(realData::defgrad_imp + d) = p.rdata(realData::deformation_gradient + d);
                    for(int d=0; d<NCOMP_TENSOR; d++) p.rdata(realData::strain_imp + d) = p.rdata(realData::strain + d);
                    p.rdata(realData::volume_imp) = p.rdata(realData::volume);

                    Real p_inf=0.0;
                    amrex::Real xp[AMREX_SPACEDIM];
                    amrex::Real gradvp[AMREX_SPACEDIM][AMREX_SPACEDIM]={0.0};
                    amrex::Real strainrate[NCOMP_TENSOR];
                    amrex::Real strain[NCOMP_TENSOR];
                    amrex::Real stress[NCOMP_TENSOR];

                    xp[XDIR]=p.pos(XDIR);
                    xp[YDIR]=p.pos(YDIR);
                    xp[ZDIR]=p.pos(ZDIR);

                    auto iv = getParticleCell(p, plo, dxi, domain);

                    lmin=(order_scheme_directional[0]==1)?0:((order_scheme_directional[0]==3 or order_scheme_directional[0]==2)?(iv[XDIR]==lo[XDIR])?0:((iv[XDIR]==hi[XDIR])?-1:-1):-1000     );
                    lmax=(order_scheme_directional[0]==1)?2:((order_scheme_directional[0]==3 or order_scheme_directional[0]==2)?(iv[XDIR]==lo[XDIR])?lmin+3:((iv[XDIR]==hi[XDIR])?lmin+3:     lmin+4):-1000);
 
                    mmin=(order_scheme_directional[1]==1)?0:((order_scheme_directional[1]==3 or order_scheme_directional[1]==2)?(iv[YDIR]==lo[YDIR])?0:((iv[YDIR]==hi[YDIR])?-1:-1):-1000     );
                    mmax=(order_scheme_directional[1]==1)?2:((order_scheme_directional[1]==3 or order_scheme_directional[1]==2)?(iv[YDIR]==lo[YDIR])?mmin+3:((iv[YDIR]==hi[YDIR])?mmin+3:     mmin+4):-1000);
 
                    nmin=(order_scheme_directional[2]==1)?0:((order_scheme_directional[2]==3 or order_scheme_directional[2]==2)?(iv[ZDIR]==lo[ZDIR])?0:((iv[ZDIR]==hi[ZDIR])?-1:-1):-1000     );
                    nmax=(order_scheme_directional[2]==1)?2:((order_scheme_directional[2]==3 or order_scheme_directional[2]==2)?(iv[ZDIR]==lo[ZDIR])?nmin+3:((iv[ZDIR]==hi[ZDIR])?nmin+3:     nmin+4):-1000);

                    if(lmin==-1000 or lmax==-1000 or mmin==-1000 or mmax==-1000 or nmin==-1000 or nmax==-1000){
                      amrex::Abort("\nError. Something wrong with min/max index values in deposit onto grid");
                    }

                    for(int n=nmin;n<nmax;n++){
                      for(int m=mmin;m<mmax;m++){
                        for(int l=lmin;l<lmax;l++){
                          IntVect ivlocal(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n);
                          if(nodalbox.contains(ivlocal)){
                              amrex::Real basisval_grad[AMREX_SPACEDIM];
                              for(int d=0;d<AMREX_SPACEDIM;d++)
                              {
                                basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme_directional,periodic,lo,hi);
                              }

                              // FIXME: Swap in updated velocity
                              gradvp[XDIR][XDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,0)*basisval_grad[XDIR];
                              gradvp[XDIR][YDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,0)*basisval_grad[YDIR];
                              gradvp[XDIR][ZDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,0)*basisval_grad[ZDIR];

                              gradvp[YDIR][XDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,1)*basisval_grad[XDIR];
                              gradvp[YDIR][YDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,1)*basisval_grad[YDIR];
                              gradvp[YDIR][ZDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,1)*basisval_grad[ZDIR];

                              gradvp[ZDIR][XDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,2)*basisval_grad[XDIR];
                              gradvp[ZDIR][YDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,2)*basisval_grad[YDIR];
                              gradvp[ZDIR][ZDIR]+=v_array(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,2)*basisval_grad[ZDIR];     

                              // gradvp[XDIR][XDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELX_INDEX)*basisval_grad[XDIR];
                              // gradvp[XDIR][YDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELX_INDEX)*basisval_grad[YDIR];
                              // gradvp[XDIR][ZDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELX_INDEX)*basisval_grad[ZDIR];

                              // gradvp[YDIR][XDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELY_INDEX)*basisval_grad[XDIR];
                              // gradvp[YDIR][YDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELY_INDEX)*basisval_grad[YDIR];
                              // gradvp[YDIR][ZDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELY_INDEX)*basisval_grad[ZDIR];

                              // gradvp[ZDIR][XDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELZ_INDEX)*basisval_grad[XDIR];
                              // gradvp[ZDIR][YDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELZ_INDEX)*basisval_grad[YDIR];
                              // gradvp[ZDIR][ZDIR]+=nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,VELZ_INDEX)*basisval_grad[ZDIR];     
                         
                          }
                        } //l loop
                      }//m loop
                    }//n loop

                    // TODO: Which particle variables is it safe to overwrite as we solve the implicit system?
                    //Calculate deformation gradient tensor at time t+dt
                    get_deformation_gradient_tensor(p,realData::defgrad_imp,gradvp,dt);

                    // Use the deformation gradient to update particle volumes
                    p.rdata(realData::jacobian) = p.rdata(realData::defgrad_imp+0)*(p.rdata(realData::defgrad_imp+4)*p.rdata(realData::defgrad_imp+8)
                                                                                            -p.rdata(realData::defgrad_imp+7)*p.rdata(realData::defgrad_imp+5))
                                                 -p.rdata(realData::defgrad_imp+1)*(p.rdata(realData::defgrad_imp+3)*p.rdata(realData::defgrad_imp+8)
                                                                                            -p.rdata(realData::defgrad_imp+6)*p.rdata(realData::defgrad_imp+5))
                                                 +p.rdata(realData::defgrad_imp+2)*(p.rdata(realData::defgrad_imp+3)*p.rdata(realData::defgrad_imp+7)
                                                                                            -p.rdata(realData::defgrad_imp+6)*p.rdata(realData::defgrad_imp+4));
                    p.rdata(realData::volume_imp) = p.rdata(realData::vol_init)*p.rdata(realData::jacobian);

                    // Update the particle strain
                    int ind=0;
                    for(int d1=0;d1<AMREX_SPACEDIM;d1++)
                    {
                      for(int d2=d1;d2<AMREX_SPACEDIM;d2++)
                      {
                        p.rdata(realData::strainrate+ind)=0.5*(gradvp[d1][d2]+gradvp[d2][d1]);
                        ind++;
                      }
                    }
                    for(int d=0;d<NCOMP_TENSOR;d++)
                    {
                      p.rdata(realData::strain_imp+d) += dt*p.rdata(realData::strainrate+d);
                    }
                    //apply axial strain
                    p.rdata(realData::strain_imp+XX) += dt*applied_strainrate;
                    p.rdata(realData::strain_imp+YY) += dt*applied_strainrate;
                    p.rdata(realData::strain_imp+ZZ) += dt*applied_strainrate;

                    // Flatten tensors
                    for(int d=0;d<NCOMP_TENSOR;d++)
                    {
                        strainrate[d]=p.rdata(realData::strainrate+d);
                        strain[d]=p.rdata(realData::strain_imp+d);
                    }

                    // Use constitutive models to calculate new particle stress
                    if(p.idata(intData::constitutive_model)==0)   //Elastic solid
                    {
                        linear_elastic(strain,strainrate,stress,p.rdata(realData::E),p.rdata(realData::nu));
                    }
                    else if(p.idata(intData::constitutive_model)==1)    //Viscous fluid with approximate EoS
                    {
                        p.rdata(realData::pressure) = p.rdata(realData::Bulk_modulus)*
                        (pow(1/p.rdata(realData::jacobian),p.rdata(realData::Gama_pressure))-1.0)+p_inf;
                        Newtonian_Fluid(strainrate,stress,p.rdata(realData::Dynamic_viscosity),p.rdata(realData::pressure));
                    }

                    for(int d=0;d<NCOMP_TENSOR;d++)
                    {
                        p.rdata(realData::stress+d)=stress[d];
                    }

                }
            });
        }

        // Loop over real and neighbor particles to update grid velocity components
        // TODO: Do we need a separare MFIter loop for this? Need to make sure that the
        //        neighbor particles we grab reflect the updated quantities. Is there any
        //        explicit sychronization needed?
        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){
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
            AMREX_GPU_DEVICE (int i) noexcept{ 
                int lmin,lmax,nmin,nmax,mmin,mmax;
                ParticleType& p = pstruct[i];
                if(p.idata(intData::phase)==0)    //Compute only for standard particles and not rigid particles with phase=1
                {
                    amrex::Real xp[AMREX_SPACEDIM];

                    xp[XDIR]=p.pos(XDIR);
                    xp[YDIR]=p.pos(YDIR);
                    xp[ZDIR]=p.pos(ZDIR);

                    auto iv = getParticleCell(p, plo, dxi, domain);

                    lmin=(order_scheme_directional[0]==1)?0:((order_scheme_directional[0]==3 or order_scheme_directional[0]==2)?(iv[XDIR]==lo[XDIR])?0:((iv[XDIR]==hi[XDIR])?-1:-1):-1000     );
                    lmax=(order_scheme_directional[0]==1)?2:((order_scheme_directional[0]==3 or order_scheme_directional[0]==2)?(iv[XDIR]==lo[XDIR])?lmin+3:((iv[XDIR]==hi[XDIR])?lmin+3:     lmin+4):-1000);
 
                    mmin=(order_scheme_directional[1]==1)?0:((order_scheme_directional[1]==3 or order_scheme_directional[1]==2)?(iv[YDIR]==lo[YDIR])?0:((iv[YDIR]==hi[YDIR])?-1:-1):-1000     );
                    mmax=(order_scheme_directional[1]==1)?2:((order_scheme_directional[1]==3 or order_scheme_directional[1]==2)?(iv[YDIR]==lo[YDIR])?mmin+3:((iv[YDIR]==hi[YDIR])?mmin+3:     mmin+4):-1000);
 
                    nmin=(order_scheme_directional[2]==1)?0:((order_scheme_directional[2]==3 or order_scheme_directional[2]==2)?(iv[ZDIR]==lo[ZDIR])?0:((iv[ZDIR]==hi[ZDIR])?-1:-1):-1000     );
                    nmax=(order_scheme_directional[2]==1)?2:((order_scheme_directional[2]==3 or order_scheme_directional[2]==2)?(iv[ZDIR]==lo[ZDIR])?nmin+3:((iv[ZDIR]==hi[ZDIR])?nmin+3:     nmin+4):-1000);

                    if(lmin==-1000 or lmax==-1000 or mmin==-1000 or mmax==-1000 or nmin==-1000 or nmax==-1000)
                    {
                      amrex::Abort("\nError. Something wrong with min/max index values in deposit onto grid");
                    }

                    for(int n=nmin;n<nmax;n++){
                      for(int m=mmin;m<mmax;m++){
                        for(int l=lmin;l<lmax;l++){
                          IntVect ivlocal(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n);
                          if(nodalbox.contains(ivlocal)){
                              amrex::Real basisval_grad[AMREX_SPACEDIM];
                              amrex::Real stress_tens[AMREX_SPACEDIM*AMREX_SPACEDIM];

                              // Flatten updated particle stress tensor
                              get_tensor(p,realData::stress,stress_tens);

                              amrex::Real basisvalue=basisval(l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme_directional,periodic,lo,hi);
                              for(int d=0;d<AMREX_SPACEDIM;d++)
                              {
                                basisval_grad[d]=basisvalder(d,l,m,n,iv[XDIR],iv[YDIR],iv[ZDIR],xp,plo,dx,order_scheme_directional,periodic,lo,hi);
                              }

                              // Calculate external forces
                              amrex::Real bforce_contrib[AMREX_SPACEDIM]=
                              {   p.rdata(realData::mass)*grav[XDIR]*basisvalue,
                                  p.rdata(realData::mass)*grav[YDIR]*basisvalue,
                                  p.rdata(realData::mass)*grav[ZDIR]*basisvalue   
                              };  
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

                              // Calculate internal force using updated particle volume
                              amrex::Real intforce_contrib[AMREX_SPACEDIM]=
                              { -p.rdata(realData::volume_imp)*tensvect[XDIR],
                                -p.rdata(realData::volume_imp)*tensvect[YDIR],
                                -p.rdata(realData::volume_imp)*tensvect[ZDIR]
                              };
                              // { -p.rdata(realData::volume)*tensvect[XDIR],
                              //   -p.rdata(realData::volume)*tensvect[YDIR],
                              //   -p.rdata(realData::volume)*tensvect[ZDIR]
                              // };

                              // Update nodal forces
                              for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                              {
                                amrex::Gpu::Atomic::AddNoRet(&nodal_data_arr(iv[XDIR]+l,iv[YDIR]+m,iv[ZDIR]+n,FRCX_INDEX+dim),
                                                             bforce_contrib[dim]+intforce_contrib[dim]);
                              }
                          }
                        }// l loop
                      }// m loop
                    }// n loop
                }
            });
        }

        // Calculate nodal velocity update using forces and nodal mass
        for ( MFIter mfi(v_data); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.tilebox();
            Box nodalbox = convert(bx, {1, 1, 1});

            const Array4<Real>& v_rhs_array = v_rhs.array(mfi);
            const Array4<Real>& nodal_arr = nodaldata.array(mfi);
      
            amrex::ParallelFor(nodalbox,[=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept { 
                if(nodal_arr(i,j,k,MASS_INDEX) > mass_tolerance){
                    v_rhs_array(i,j,k,0) = nodal_arr(i,j,k,FRCX_INDEX+0)/nodal_arr(i,j,k,MASS_INDEX);
                    v_rhs_array(i,j,k,1) = nodal_arr(i,j,k,FRCX_INDEX+1)/nodal_arr(i,j,k,MASS_INDEX);
                    v_rhs_array(i,j,k,2) = nodal_arr(i,j,k,FRCX_INDEX+2)/nodal_arr(i,j,k,MASS_INDEX);
                }
                else{
                    v_rhs_array(i,j,k,0) = 0.0;
                    v_rhs_array(i,j,k,1) = 0.0;
                    v_rhs_array(i,j,k,2) = 0.0;
                }
            });
        }
    };  
    /* ----------------------------------------------------------------*
     *                        END OF RHS KERNEL                        *
     * ----------------------------------------------------------------*/

    // Define the post-update function
    auto post_update_function = [&](Vector<MultiFab>& S_data, const Real /* time */) {
        // fill periodic ghost cells
        // S_data[0].FillBoundary(geom.periodicity());
    };

    // Set up the AMReX time integrator
    TimeIntegrator<Vector<MultiFab> > integrator(v_vect_old);
    integrator.set_rhs(rhs_function);
    integrator.set_post_update(post_update_function);
    
    // Advance from state_old at time to state_new at time + dt
    integrator.advance(v_vect_old, v_vect_new, time, dt);    

    // Copy updated velocity back into nodal data MF
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        Box nodalbox = convert(box, {1, 1, 1});
 
        Array4<Real> nodal_data_arr=nodaldata.array(mfi);
        Array4<Real> vnew_data_arr = v_vect_new[0].array(mfi);
        Array4<Real> vold_data_arr = v_vect_old[0].array(mfi);

        amrex::ParallelFor(nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // FIXME: Do we still need to account for small nodal masses when using implicit solve?
            if(nodal_data_arr(i,j,k,MASS_INDEX)>=mass_tolerance)
            {
                for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                {
                    nodal_data_arr(i,j,k,VELX_INDEX + dim) = vnew_data_arr(i,j,k,dim);
                }
            }
            else
            {
                for(int dim=0;dim<AMREX_SPACEDIM;dim++)
                {
                    nodal_data_arr(i,j,k,VELX_INDEX + dim) = 0.0;
                }
            }
        });
    } 
}
