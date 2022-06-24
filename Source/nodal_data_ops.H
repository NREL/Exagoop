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
void initialise_shape_function_indices(amrex::iMultiFab &shapefunctionindex,const amrex::Geometry geom);
void nodal_bcs(const amrex::Geometry geom,amrex::MultiFab &nodaldata,const amrex::Real& dt);

#endif
