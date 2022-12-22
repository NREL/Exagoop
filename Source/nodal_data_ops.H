#ifndef NODAL_DATA_OPS_H_
#define NODAL_DATA_OPS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <constants.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

void write_plot_file(std::string fname, amrex::MultiFab &nodaldata, amrex::Vector<std::string> fieldnames, 
                           amrex::Geometry geom, amrex::BoxArray ba, amrex::DistributionMapping dm,amrex::Real time);

void backup_current_velocity(amrex::MultiFab &nodaldata);
void store_delta_velocity(amrex::MultiFab &nodaldata);
void nodal_update(amrex::MultiFab &nodaldata,const amrex::Real& dt, const amrex::Real& mass_tolerance);
void nodal_detect_contact(amrex::MultiFab &nodaldata,const amrex::Geometry geom,amrex::Real& contact_tolerance,amrex::Real vely);
void initialise_shape_function_indices(amrex::iMultiFab &shapefunctionindex,const amrex::Geometry geom);


void nodal_bcs(const amrex::Geometry geom,amrex::MultiFab &nodaldata,
int bclo[AMREX_SPACEDIM],int bchi[AMREX_SPACEDIM],
amrex::Real wall_mu_lo[AMREX_SPACEDIM],amrex::Real wall_mu_hi[AMREX_SPACEDIM],
amrex::Real wall_vel_lo[AMREX_SPACEDIM*AMREX_SPACEDIM],
amrex::Real wall_vel_hi[AMREX_SPACEDIM*AMREX_SPACEDIM],
const amrex::Real& dt);
void nodal_bcs(const amrex::Geometry geom, amrex::MultiFab &nodaldata,const amrex::Real& dt);
void CalculateSurfaceIntegralOnBG(const amrex::Geometry geom,amrex::MultiFab &nodaldata, int nodaldataindex,amrex::Real &integral_value);
void CalculateInterpolationError(const amrex::Geometry geom,amrex::MultiFab &nodaldata, int nodaldataindex);

void nodal_levelset_bcs(amrex::MultiFab &nodaldata,const amrex::Geometry geom,
        amrex::Real& dt,int lsetbc,amrex::Real lset_wall_mu);

#endif
