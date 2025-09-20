// clang-format off
#include <nodal_data_ops.H>
#include <mpm_eb.H>
#include <mpm_kernels.H>
// clang-format on

using namespace amrex;

void write_grid_file(std::string fname, MultiFab &nodaldata,
                     Vector<std::string> fieldnames, Geometry geom, BoxArray ba,
                     DistributionMapping dm, Real time) {
  MultiFab plotmf(ba, dm, nodaldata.nComp(), 0);
  average_node_to_cellcenter(plotmf, 0, nodaldata, 0, nodaldata.nComp());
  WriteSingleLevelPlotfile(fname, plotmf, fieldnames, geom, time, 0);
}

void backup_current_velocity(MultiFab &nodaldata) {
  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<Real> nodal_data_arr = nodaldata.array(mfi);

    amrex::ParallelFor(nodalbox,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         if (nodal_data_arr(i, j, k, MASS_INDEX) > zero) {
                           nodal_data_arr(i, j, k, MASS_OLD_INDEX) =
                               nodal_data_arr(i, j, k, MASS_INDEX);
                           for (int d = 0; d < AMREX_SPACEDIM; d++) {
                             // Storing V(t) values at these memory locations.
                             // This is to be followed by a call to
                             // store_delta_velocity
                             nodal_data_arr(i, j, k, DELTA_VELX_INDEX + d) =
                                 nodal_data_arr(i, j, k, VELX_INDEX + d);
                           }
                         }
                       });
  }
}

void nodal_levelset_bcs(MultiFab &nodaldata, const Geometry geom,
                        amrex::Real &dt, int lsetbc, amrex::Real lset_wall_mu) {
  // need something more sophisticated
  // but lets get it working!
  //
  int lsref = mpm_ebtools::ls_refinement;
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const auto dx = geom.CellSizeArray();

  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<Real> nodal_data_arr = nodaldata.array(mfi);
    Array4<Real> lsarr = mpm_ebtools::lsphi->array(mfi);

    amrex::ParallelFor(
        nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          IntVect nodeid(i, j, k);
          IntVect refined_nodeid(i * lsref, j * lsref, k * lsref);

          if (lsarr(refined_nodeid) < TINYVAL &&
              nodal_data_arr(nodeid, MASS_INDEX) > zero) {
            Real relvel_in[AMREX_SPACEDIM];
            Real relvel_out[AMREX_SPACEDIM];
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              relvel_in[d] = nodal_data_arr(nodeid, VELX_INDEX + d);
              relvel_out[d] = relvel_in[d];
            }

            amrex::Real eps = 0.00001;
            amrex::Real xp[AMREX_SPACEDIM] = {plo[XDIR] + i * dx[XDIR],
                                              plo[YDIR] + j * dx[YDIR],
                                              plo[ZDIR] + k * dx[ZDIR]};

            // amrex::Real dist=get_levelset_value(lsarr,plo,dx,xp,lsref);
            // dist=amrex::Math::abs(dist);

            amrex::Real normaldir[AMREX_SPACEDIM] = {1.0, 0.0, 0.0};

            get_levelset_grad(lsarr, plo, dx, xp, lsref, normaldir);
            amrex::Real gradmag = std::sqrt(normaldir[XDIR] * normaldir[XDIR] +
                                            normaldir[YDIR] * normaldir[YDIR] +
                                            normaldir[ZDIR] * normaldir[ZDIR]);

            normaldir[XDIR] = normaldir[XDIR] / (gradmag + TINYVAL);
            normaldir[YDIR] = normaldir[YDIR] / (gradmag + TINYVAL);
            normaldir[ZDIR] = normaldir[ZDIR] / (gradmag + TINYVAL);

            int modify_pos =
                applybc(relvel_in, relvel_out, lset_wall_mu, normaldir, lsetbc);

            nodal_data_arr(nodeid, VELX_INDEX + XDIR) = relvel_out[XDIR];
            nodal_data_arr(nodeid, VELX_INDEX + YDIR) = relvel_out[YDIR];
            nodal_data_arr(nodeid, VELX_INDEX + ZDIR) = relvel_out[ZDIR];
          }
        });
  }
}

void store_delta_velocity(MultiFab &nodaldata) {
  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<Real> nodal_data_arr = nodaldata.array(mfi);

    amrex::ParallelFor(nodalbox,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         if (nodal_data_arr(i, j, k, MASS_INDEX) > zero) {
                           for (int d = 0; d < AMREX_SPACEDIM; d++) {
                             nodal_data_arr(i, j, k, DELTA_VELX_INDEX + d) =
                                 nodal_data_arr(i, j, k, VELX_INDEX + d) -
                                 nodal_data_arr(i, j, k, DELTA_VELX_INDEX + d);
                           }
                         }
                       });
  }
}

void nodal_update(MultiFab &nodaldata, const amrex::Real &dt,
                  const amrex::Real &mass_tolerance) {
  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<Real> nodal_data_arr = nodaldata.array(mfi);

    amrex::ParallelFor(
        nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (nodal_data_arr(i, j, k, MASS_INDEX) >= mass_tolerance) {
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              nodal_data_arr(i, j, k, VELX_INDEX + d) +=
                  nodal_data_arr(i, j, k, FRCX_INDEX + d) /
                  nodal_data_arr(i, j, k, MASS_INDEX) * dt;
            }
          } else {
            nodal_data_arr(i, j, k, VELX_INDEX) = 0.0;
            nodal_data_arr(i, j, k, VELY_INDEX) = 0.0;
            nodal_data_arr(i, j, k, VELZ_INDEX) = 0.0;
          }
        });
  }
}

void nodal_detect_contact(
    MultiFab &nodaldata, const Geometry geom, amrex::Real &contact_tolerance,
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>,
                    numrigidbodies>
        velocity) {

  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<Real> nodal_data_arr = nodaldata.array(mfi);

    amrex::ParallelFor(
        nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (nodal_data_arr(i, j, k, MASS_INDEX) > contact_tolerance and
              nodal_data_arr(i, j, k, MASS_RIGID_INDEX) > contact_tolerance and
              int(nodal_data_arr(i, j, k, RIGID_BODY_ID)) != -1) {
            amrex::Real contact_alpha = 0.0;
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              contact_alpha +=
                  (nodal_data_arr(i, j, k, VELX_INDEX + d) -
                   velocity[int(nodal_data_arr(i, j, k, RIGID_BODY_ID))][d]) *
                  nodal_data_arr(i, j, k, NORMALX + d);
            }

            if (contact_alpha >= 0) {
              amrex::Real V_relative = 0.0;
              for (int d = 0; d < AMREX_SPACEDIM; d++) {
                V_relative +=
                    (nodal_data_arr(i, j, k, VELX_INDEX + d) -
                     velocity[int(nodal_data_arr(i, j, k, RIGID_BODY_ID))][d]) *
                    nodal_data_arr(i, j, k, NORMALX + d);
              }
              for (int d = 0; d < AMREX_SPACEDIM; d++) {
                nodal_data_arr(i, j, k, VELX_INDEX + d) -=
                    V_relative * (nodal_data_arr(i, j, k, NORMALX + d));
              }
            } else {
              // nodal_data_arr(i,j,k,VELX_INDEX+1)=0.0;
            }
          }
        });
  }
}

void initialise_shape_function_indices(iMultiFab &shapefunctionindex,
                                       const amrex::Geometry geom) {
  const int *domloarr = geom.Domain().loVect();
  const int *domhiarr = geom.Domain().hiVect();

  int periodic[AMREX_SPACEDIM] = {geom.isPeriodic(XDIR), geom.isPeriodic(YDIR),
                                  geom.isPeriodic(ZDIR)};

  GpuArray<int, AMREX_SPACEDIM> lo = {domloarr[0], domloarr[1], domloarr[2]};
  GpuArray<int, AMREX_SPACEDIM> hi = {domhiarr[0], domhiarr[1], domhiarr[2]};

  for (MFIter mfi(shapefunctionindex); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<int> shapefunctionindex_arr = shapefunctionindex.array(mfi);

    amrex::ParallelFor(nodalbox,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         if (i == lo[0]) {
                           shapefunctionindex_arr(i, j, k, 0) = 0;
                         } else if (i == lo[0] + 1) {
                           shapefunctionindex_arr(i, j, k, 0) = 1;
                         } else if (i == hi[0]) {
                           shapefunctionindex_arr(i, j, k, 0) = 4;
                         } else if (i == hi[0] - 1) {
                           shapefunctionindex_arr(i, j, k, 0) = 3;
                         } else {
                           shapefunctionindex_arr(i, j, k, 0) = 2;
                         }

                         if (j == lo[1]) {
                           shapefunctionindex_arr(i, j, k, 1) = 0;
                         } else if (j == lo[1] + 1) {
                           shapefunctionindex_arr(i, j, k, 1) = 1;
                         } else if (j == hi[1]) {
                           shapefunctionindex_arr(i, j, k, 1) = 4;
                         } else if (j == hi[1] - 1) {
                           shapefunctionindex_arr(i, j, k, 1) = 3;
                         } else {
                           shapefunctionindex_arr(i, j, k, 1) = 2;
                         }

                         if (k == lo[2]) {
                           shapefunctionindex_arr(i, j, k, 2) = 0;
                         } else if (k == lo[2] + 1) {
                           shapefunctionindex_arr(i, j, k, 2) = 1;
                         } else if (k == hi[2]) {
                           shapefunctionindex_arr(i, j, k, 2) = 4;
                         } else if (k == hi[2] - 1) {
                           shapefunctionindex_arr(i, j, k, 2) = 3;
                         } else {
                           shapefunctionindex_arr(i, j, k, 2) = 2;
                         }
                       });
  }
}

void nodal_bcs(const amrex::Geometry geom, MultiFab &nodaldata,
               int bcloarr[AMREX_SPACEDIM], int bchiarr[AMREX_SPACEDIM],
               Real wall_mu_loarr[AMREX_SPACEDIM],
               Real wall_mu_hiarr[AMREX_SPACEDIM],
               Real wall_vel_loarr[AMREX_SPACEDIM * AMREX_SPACEDIM],
               Real wall_vel_hiarr[AMREX_SPACEDIM * AMREX_SPACEDIM],
               const amrex::Real &dt) {

  const int *domloarr = geom.Domain().loVect();
  const int *domhiarr = geom.Domain().hiVect();

  int periodic[AMREX_SPACEDIM] = {geom.isPeriodic(XDIR), geom.isPeriodic(YDIR),
                                  geom.isPeriodic(ZDIR)};

  // gpu capture stuff
  GpuArray<int, AMREX_SPACEDIM> domlo = {domloarr[0], domloarr[1], domloarr[2]};
  GpuArray<int, AMREX_SPACEDIM> domhi = {domhiarr[0], domhiarr[1], domhiarr[2]};

  int bclo[] = {bcloarr[0], bcloarr[1], bcloarr[2]};
  int bchi[] = {bchiarr[0], bchiarr[1], bchiarr[2]};

  Real wall_mu_lo[] = {wall_mu_loarr[0], wall_mu_loarr[1], wall_mu_loarr[2]};
  Real wall_mu_hi[] = {wall_mu_hiarr[0], wall_mu_hiarr[1], wall_mu_hiarr[2]};

  GpuArray<Real, AMREX_SPACEDIM * AMREX_SPACEDIM> wall_vel_lo;
  GpuArray<Real, AMREX_SPACEDIM * AMREX_SPACEDIM> wall_vel_hi;

  for (int d = 0; d < AMREX_SPACEDIM * AMREX_SPACEDIM; d++) {
    wall_vel_lo[d] = wall_vel_loarr[d];
    wall_vel_hi[d] = wall_vel_hiarr[d];
  }

  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<Real> nodal_data_arr = nodaldata.array(mfi);

    amrex::ParallelFor(
        nodalbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          IntVect nodeid(i, j, k);
          Real relvel_in[AMREX_SPACEDIM], relvel_out[AMREX_SPACEDIM];
          Real wallvel[AMREX_SPACEDIM] = {0.0, 0.0, 0.0};

          for (int d = 0; d < AMREX_SPACEDIM; d++) {
            relvel_in[d] = nodal_data_arr(nodeid, VELX_INDEX + d);
            relvel_out[d] = relvel_in[d];
          }

          if (nodeid[XDIR] == domlo[XDIR]) {
            int dir = XDIR;
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              wallvel[d] = wall_vel_lo[dir * AMREX_SPACEDIM + d];
              relvel_in[d] -= wallvel[d];
            }

            Real normaldir[AMREX_SPACEDIM] = {1.0, 0.0, 0.0};
            int tmp = applybc(relvel_in, relvel_out, wall_mu_lo[XDIR],
                              normaldir, bclo[XDIR]);
            // amrex::Print()<<"\nVelx 0 = "<<relvel_out[0];*/

          } else if (nodeid[XDIR] == (domhi[XDIR] + 1)) {
            int dir = XDIR;
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              wallvel[d] = wall_vel_hi[dir * AMREX_SPACEDIM + d];
              relvel_in[d] -= wallvel[d];
            }

            Real normaldir[AMREX_SPACEDIM] = {-1.0, 0.0, 0.0};
            int tmp = applybc(relvel_in, relvel_out, wall_mu_hi[XDIR],
                              normaldir, bchi[XDIR]);
            // relvel_out[0]=0.0;
            // amrex::Print()<<"\nVelx 1 = "<<relvel_out[0];
          } else if (nodeid[YDIR] == domlo[YDIR]) {
            int dir = YDIR;
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              wallvel[d] = wall_vel_lo[dir * AMREX_SPACEDIM + d];
              relvel_in[d] -= wallvel[d];
            }

            Real normaldir[AMREX_SPACEDIM] = {0.0, 1.0, 0.0};
            int tmp = applybc(relvel_in, relvel_out, wall_mu_lo[YDIR],
                              normaldir, bclo[YDIR]);
            // relvel_out[1]=0.0;
          } else if (nodeid[YDIR] == (domhi[YDIR] + 1)) {
            int dir = YDIR;
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              wallvel[d] = wall_vel_hi[dir * AMREX_SPACEDIM + d];
              relvel_in[d] -= wallvel[d];
            }

            Real normaldir[AMREX_SPACEDIM] = {0.0, -1.0, 0.0};
            int tmp = applybc(relvel_in, relvel_out, wall_mu_hi[YDIR],
                              normaldir, bchi[YDIR]);
            // relvel_out[1]=0.0;
          } else if (nodeid[ZDIR] == domlo[ZDIR]) {
            int dir = ZDIR;
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              wallvel[d] = wall_vel_lo[dir * AMREX_SPACEDIM + d];
              relvel_in[d] -= wallvel[d];
            }

            Real normaldir[AMREX_SPACEDIM] = {0.0, 0.0, 1.0};
            int tmp = applybc(relvel_in, relvel_out, wall_mu_lo[ZDIR],
                              normaldir, bclo[ZDIR]);
            // relvel_out[2]=0.0;
          } else if (nodeid[ZDIR] == (domhi[ZDIR] + 1)) {
            int dir = ZDIR;
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
              wallvel[d] = wall_vel_hi[dir * AMREX_SPACEDIM + d];
              relvel_in[d] -= wallvel[d];
            }

            Real normaldir[AMREX_SPACEDIM] = {0.0, 0.0, -1.0};
            int tmp = applybc(relvel_in, relvel_out, wall_mu_hi[ZDIR],
                              normaldir, bchi[ZDIR]);
            // relvel_out[2]=0.0;
          } else // nothing to do
          {
          }

          for (int d = 0; d < AMREX_SPACEDIM; d++) {
            nodal_data_arr(nodeid, VELX_INDEX + d) = relvel_out[d] + wallvel[d];
          }
          // amrex::Print()<<"\nX = Vel "<<nodal_data_arr(nodeid,0);
        });
  }
}

void CalculateSurfaceIntegralOnBG(const amrex::Geometry geom,
                                  amrex::MultiFab &nodaldata,
                                  int nodaldataindex,
                                  amrex::Real &integral_value) {
  // For the time being, calculate the surface integral of y-velocity on the
  // bottom surface ie y=0 and validate the code
  integral_value = 0.0;

  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<Real> nodal_data_arr = nodaldata.array(mfi);

    // amrex::ParallelFor(nodalbox,[=,&integral_value]AMREX_GPU_DEVICE (int
    // i,int j,int k) noexcept
    //{
    //	integral_value+=(nodal_data_arr(i,j,k,nodaldataindex)+nodal_data_arr(i+1,j,k,nodaldataindex)+nodal_data_arr(i,j,k+1,nodaldataindex)+nodal_data_arr(i+1,j,k+1,nodaldataindex))/4.0*dx[0]*dx[2];
    //	});
  }
#ifdef BL_USE_MPI
  ParallelDescriptor::ReduceRealSum(integral_value);
#endif
}

void CalculateInterpolationError(const amrex::Geometry geom,
                                 amrex::MultiFab &nodaldata,
                                 int nodaldataindex) {
  const int *domloarr = geom.Domain().loVect();
  const int *domhiarr = geom.Domain().hiVect();
  const auto dx = geom.CellSizeArray();
  const auto Pi = 4.0 * atan(1.0);

  for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    Box nodalbox = convert(bx, {1, 1, 1});

    Array4<Real> nodal_data_arr = nodaldata.array(mfi);

    amrex::ParallelFor(nodalbox,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         nodal_data_arr(i, j, k, nodaldataindex) =
                             nodal_data_arr(i, j, k, VELX_INDEX) -
                             cos(2.0 * Pi * (domloarr[0] + i * dx[0]) / 1.0);
                       });
  }
}
