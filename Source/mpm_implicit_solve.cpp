#include <mpm_particle_container.H>
#include <mpm_implicit.H>

void MPMParticleContainer::implicitUpdate(int implicit_solve, 
                                          Vector<MultiFab>& v_vect_old, 
                                          Vector<MultiFab>& v_vect_new, 
                                          MultiFab& fext, 
                                          MultiFab& nodaldata,
                                          Real mass_tolerance,
                                          Real time, Real dt){

    if(implicit_solve == 1){
        implicitSolveM1(v_vect_old, v_vect_new, fext, nodaldata, mass_tolerance, time, dt);    
    } else if(implicit_solve == 2){
    
    } else {
        Error("ERROR: Invalid implicit method selected!\n");
    }
}

void MPMParticleContainer::implicitSolveM1(Vector<MultiFab>& v_vect_old, 
                                           Vector<MultiFab>& v_vect_new,
                                           MultiFab& fext, 
                                           MultiFab& nodaldata, 
                                           Real mass_tolerance, 
                                           Real time, Real dt){

    const int lev = 0;

    // Define RHS function for the AMReX integrator
    auto rhs_function = [&](Vector<MultiFab>& S_rhs,
                      const Vector<MultiFab>& S_data, const Real /* time */) {

        auto& v_data = S_data[0];
        auto& v_rhs  = S_rhs[0];

        // Calculate internal forces
  
        // Add on external forces
        for ( MFIter mfi(v_data); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.tilebox();
            Box nodalbox = convert(bx, {1, 1, 1});

            const Array4<const Real>& v_array = v_data.array(mfi);
            const Array4<Real>& v_rhs_array = v_rhs.array(mfi);
            const Array4<Real>& fe_array = fext.array(mfi);
      
            amrex::ParallelFor(nodalbox,[=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept { 
                v_rhs_array(i,j,k,0) = fe_array(i,j,k,0);
                v_rhs_array(i,j,k,1) = fe_array(i,j,k,1);
                v_rhs_array(i,j,k,2) = fe_array(i,j,k,2);
            });
        }
    };

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
            for(int dim=0;dim<AMREX_SPACEDIM;dim++)
            {
                if(nodal_data_arr(i,j,k,MASS_INDEX)>=mass_tolerance)
                {
                    vnew_data_arr(i,j,k,dim)/=nodal_data_arr(i,j,k,MASS_INDEX);
                }
                else
                {
                    vnew_data_arr(i,j,k,dim) = 0.0;
                }
                // amrex::Real tempV = (nodal_data_arr(i,j,k,MASS_INDEX) > mass_tolerance) ? vold_data_arr(i,j,k,dim) + nodal_data_arr(i,j,k,FRCX_INDEX+dim)/nodal_data_arr(i,j,k,MASS_INDEX)*dt:0;
                // nodal_data_arr(i,j,k,VELX_INDEX+dim) = vnew_data_arr(i,j,k,dim);
                // printf("nodal data(%i, %i, %i, %i) = %.6e, vs %.6e\n", i, j, k, dim, vnew_data_arr(i,j,k,dim), tempV);
            }
        });
    } 
}


// intheck_flag(void* flagvalue, const char* funcname, int opt){
//   if (opt == 0 && flagvalue == nullptr) {
//     if (amrex::ParallelDescriptor::IOProcessor()) {
//       fprintf(
//         stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
//         funcname);
//       amrex::Abort("abort");
//     }
//     return (1);
//   }
//   if (opt == 1) {
//     int* errflag = static_cast<int*>(flagvalue);
//     if (*errflag < 0) {
//       if (amrex::ParallelDescriptor::IOProcessor()) {
//         fprintf(
//           stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname,
//           *errflag);
//         amrex::Abort("abort");
//       }
//       return (1);
//     }
//   } else if (opt == 2 && flagvalue == nullptr) {
//     if (amrex::ParallelDescriptor::IOProcessor()) {
//       fprintf(
//         stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
//         funcname);
//       amrex::Abort("abort");
//     }
//     return (1);
//   }
// 
//   return (0);
// }
// 
// int CvodeIntegrator::init(int ncell_x, int ncell_y, int ncell_z){
// 
//     // Parse "ode" input parameters
//     amrex::ParmParse pp("ode");
//     pp.query("verbose", verbose);
//     pp.query("rtol", relTol);
//     pp.query("atol", absTol);
//     pp.query("max_nls_iters", max_nls_iters);
// 
//     /* --------------------------------------
//        Initialize cvode objects and user data
//        -------------------------------------- */
// 
//     // Allocate solution vector
//     npx = ncell_x; npy = ncell_y; npz = ncell_z;
//     np = npx * npy * npz;
//     int neq_tot = 3 * np;
//     y = N_VNew_Serial(neq_tot, sunctx);
//     if (check_flag(static_cast<void*>(y), "N_VNew_Serial", 0) != 0) {
//         return (1);
//     }
// 
//     // Using BDF method for time discretization
//     cvode_mem = CVodeCreate(CV_BDF, sunctx);
//     if (check_flag(static_cast<void*>(cvode_mem), "CVodeCreate", 0) != 0) {
//         return (1);
//     }
// 
//     // Create and allocate user data
//     udata_g = new CVODEUserData{};
//     allocUserData(udata_g);
//     if (check_flag(static_cast<void*>(udata_g), "allocUserData", 2) != 0) {
//         return (1);
//     }
//     int flag = CVodeSetUserData(cvode_mem, udata_g);
//     if (check_flag(&flag, "CVodeSetUserData", 1) != 0) {
//         return (1);
//     }
// 
//     // Set initial solver time, soln vector and RHS function
//     amrex::Real time = 0.0;
//     flag = CVodeInit(cvode_mem, cF_RHS, time, y);
//     if (check_flag(&flag, "CVodeInit", 1) != 0) {
//         return (1);
//     }
// 
//     // Setting scalar tolerances for now, may be an issue for multiphase flows, etc.  
//     flag = CVodeSStolerances(cvode_mem, relTol, absTol);
//     if (check_flag(&flag, "CVodeInit", 1) != 0) {
//         return (1);
//     }
// 
//     // For now just use GMRES
//     LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
//     if (check_flag(static_cast<void*>(LS), "SUNLinSol_SPGMR", 0) != 0) {
//         return (1);
//     }
//     flag = CVodeSetLinearSolver(cvode_mem, LS, nullptr);
//     if (check_flag(&flag, "CVodeSetLinearSolver", 1) != 0) {
//       return (1);
//     }
// 
//     // Using mostly CVODE deffault runtime options
//     flag = CVodeSetMaxNonlinIters(cvode_mem, max_nls_iters); // Max newton iter.
//     if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1) != 0) {
//         return (1);
//     }
//     // flag = CVodeSetMaxErrTestFails(cvode_mem, 100); // Max Err.test failure
//     // if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1) != 0) {
//     //     return (1);
//     // }
//     // flag = CVodeSetErrHandlerFn(
//     //   cvode_mem, cvode::cvodeErrHandler, nullptr); // Err. handler funct.
//     // if (utils::check_flag(&flag, "CVodeSetErrHandlerFn", 1) != 0) {
//     //     return (1);
//     // }
//     // flag = CVodeSetMaxNumSteps(cvode_mem, 10000); // Max substeps
//     // if (check_flag(&flag, "CVodeSetMaxNumSteps", 1) != 0) {
//     //     return (1);
//     // }
//     // flag = CVodeSetMaxOrd(cvode_mem, udata_g->maxOrder); // Max order
//     // if (utils::check_flag(&flag, "CVodeSetMaxOrd", 1) != 0) {
//     //     return (1);
//     // }
// 
//     return(0);
// }
// 
// void CvodeIntegrator::allocUserData(CVODEUserData* udata){
// 
//     udata->npts = np;
//     udata->npts_x = npx;
//     udata->npts_y = npy;
//     udata->npts_z = npz;
//     udata->fext = new amrex::Real[3*np];
// }
// 
// void CvodeIntegrator::solve(amrex::Real ctime, amrex::Real dt, amrex::MultiFab& nodaldata, amrex::MultiFab& f_ext, amrex::Real mass_tolerance){
// 
//     // NOTE: assuming 3D case for now...
// 
//     // Get pointer to soln. data
//     amrex::Real* yvec_d = N_VGetArrayPointer(y);
// 
//     // Fill soln. vector with velocity components at t=n
//     for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
//     {
//         const Box& bx=mfi.validbox();
//         Box nodalbox = convert(bx, {1, 1, 1});
// 
//         Array4<Real> nodal_data_arr=nodaldata.array(mfi);
//         Array4<Real> fext_data_arr=f_ext.array(mfi);
// 
//         amrex::ParallelFor(nodalbox,[=]
//         AMREX_GPU_DEVICE (int i,int j,int k) noexcept
//         {
//             int idx = 3*(k*(npx*npy) + j*(npx) + i);
// 
//             // Flatten grid velocity data
//             yvec_d[idx + 0] = nodal_data_arr(i,j,k,VELX_INDEX + 0);
//             yvec_d[idx + 1] = nodal_data_arr(i,j,k,VELX_INDEX + 1);
//             yvec_d[idx + 2] = nodal_data_arr(i,j,k,VELX_INDEX + 2);
//     
//             // Flatten grid external forces
//             // FIXME: Do we still need to worry about small mass values w/ implicit system?
//             if(nodal_data_arr(i,j,k,MASS_INDEX) >= mass_tolerance){
//                 udata_g->fext[idx + 0] = fext_data_arr(i,j,k,0) / nodal_data_arr(i,j,k,MASS_INDEX);
//                 udata_g->fext[idx + 1] = fext_data_arr(i,j,k,1) / nodal_data_arr(i,j,k,MASS_INDEX);
//                 udata_g->fext[idx + 2] = fext_data_arr(i,j,k,2) / nodal_data_arr(i,j,k,MASS_INDEX);
//             } else { 
//                 udata_g->fext[idx + 0] = 0.0;
//                 udata_g->fext[idx + 1] = 0.0;
//                 udata_g->fext[idx + 2] = 0.0;
//             }
//         });
//     }
// 
//     // Reinitialize cvode with current time
//     CVodeReInit(cvode_mem, ctime, y);
// 
//     // Solve NL system
//     CVode(cvode_mem, ctime + dt, y, &CvodeActual_time_final, CV_NORMAL);
// 
//     if(verbose > 0){
//         long int nfe = 0;
//         long int nfeLS = 0;
//         CVodeGetNumRhsEvals(cvode_mem, &nfe);
//         if (LS != nullptr) {
//           CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
//         }
//         const long int nfe_tot = nfe + nfeLS;
// 
//         long int nnli = 0;
//         long int nli = 0;
//         CVodeGetNumNonlinSolvIters(cvode_mem, &nnli);
//         if(LS != nullptr){
//             CVodeGetNumLinIters(cvode_mem, &nli);
//         }
//         amrex::Print() << "\nTotal RHS fcn evals: " << nfe_tot << ",\t\t Nonlinear iterations: " << nnli << ",\t\t Linear iteration: " << nli;
//     }
// 
//     // Add updated velocities back into nodal MF
//     for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
//     {
//         const Box& bx=mfi.validbox();
//         Box nodalbox = convert(bx, {1, 1, 1});
// 
//         Array4<Real> nodal_data_arr=nodaldata.array(mfi);
// 
//         amrex::ParallelFor(nodalbox,[=]
//         AMREX_GPU_DEVICE (int i,int j,int k) noexcept
//         {
//             int idx = 3*(k*(npx*npy) + j*(npx) + i);
// 
//             // Unflatten grid velocity data
//             // FIXME: Do we still need to worry about small mass values w/ implicit system?
//             if(nodal_data_arr(i,j,k,MASS_INDEX) >= mass_tolerance){
//                 nodal_data_arr(i,j,k,VELX_INDEX + 0) = yvec_d[idx + 0];
//                 nodal_data_arr(i,j,k,VELX_INDEX + 1) = yvec_d[idx + 1];
//                 nodal_data_arr(i,j,k,VELX_INDEX + 2) = yvec_d[idx + 2];
//             } else {
//                 nodal_data_arr(i,j,k,VELX_INDEX + 0) = 0.0;
//                 nodal_data_arr(i,j,k,VELX_INDEX + 1) = 0.0;
//                 nodal_data_arr(i,j,k,VELX_INDEX + 2) = 0.0;
//             }
//   
//         });
//     }
// }
// 
// int CvodeIntegrator::cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, void* user_data){
// 
//     // Get pointers to y and ydot 
//     amrex::Real* yvec_d = N_VGetArrayPointer(y_in);
//     amrex::Real* ydot_d = N_VGetArrayPointer(ydot_in);
// 
//     // Unpack user data
//     auto* udata = static_cast<CVODEUserData*>(user_data);    
//     int ncells = udata->npts;
//     int nx = udata->npts_x;
//     int ny = udata->npts_y;
//     int nz = udata->npts_z;
//     auto* fe = udata->fext;
// 
//     // Domain data
//     const int lev = 0;
// 
//     // Loop over particles to get internal force contribution for each grid node
//     // FIXME: Is it ok to have MFIter in CVODE kernel?
//     for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
//     {
//         const amrex::Box& box = mfi.tilebox();
//         Box nodalbox = convert(box, {1, 1, 1});
// 
//         int gid = mfi.index();
//         int tid = mfi.LocalTileIndex();
//         auto index = std::make_pair(gid, tid);
// 
//         // auto& ptile = plev[index];
//         // auto& aos   = ptile.GetArrayOfStructs();
//         // int np = aos.numRealParticles();
//         // int ng =aos.numNeighborParticles();
//         // int nt = np+ng;
// 
//         // Array4<Real> nodal_data_arr=nodaldata.array(mfi);
//         // Array4<Real> fext_data_arr = f_ext.array(mfi);
// 
//         // ParticleType* pstruct = aos().dataPtr();
//     }
// 
// 
//     // Add on external force contribution
//     amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept {
//         ydot_d[icell*3 + 0] = fe[icell*3 + 0];
//         ydot_d[icell*3 + 1] = fe[icell*3 + 1];
//         ydot_d[icell*3 + 2] = fe[icell*3 + 2];
//     });
// 
//     return(0);
// }




