#include <AMReX_TimeIntegrator.H>

// struct CVODEUserData
// {
// 
//   int npts;                         // Number of points
//   int npts_x;                       // Number of points in x-dir
//   int npts_y;                       // Number of points in y-dir
//   int npts_z;                       // Number of points in z-dir
//   amrex::Real* fext = nullptr;      // External force array
// 
// };
// 
// class CvodeIntegrator
// {
// public:
// 
//     int init(int ncell_x, int ncell_y, int ncell_z);
//     void solve(amrex::Real ctime, amrex::Real dt, amrex::MultiFab& nodaldata, amrex::MultiFab& f_ext, amrex::Real mass_tolerance);
//     static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, void* user_data);
// 
// private: 
//       
//     // User input variables
//     int verbose = 0;
//     amrex::Real relTol = 1.0e-6;
//     amrex::Real absTol = 1.0e-6;
//     int max_nls_iters = 5;        
// 
//     // Problem size variables
//     int np;
//     int npx;
//     int npy;
//     int npz;
// 
//     // Sundials objects
//     SUNContext sunctx;
//     void* cvode_mem = nullptr;
//     CVODEUserData* udata_g = nullptr;
//     SUNLinearSolver LS = nullptr;
// 
//     // Sundials vectors
//     N_Vector y = nullptr;
// 
//     // Misc. Sundials variables
//     amrex::Real CvodeActual_time_final;
// 
//     void allocUserData(CVODEUserData* udata);
// 
// };
