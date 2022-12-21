#ifndef _MPM_EB_H_
#define _MPM_EB_H_

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2_IF_Cylinder.H>
<<<<<<< HEAD
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Translation.H>

#include <AMReX_EB_utils.H>
#include <constants.H>
#include <interpolants.H>

using namespace amrex;

namespace mpm_ebtools
{
    extern EBFArrayBoxFactory* ebfactory;
    extern MultiFab *lsphi;
    extern int ls_refinement;
    extern bool using_levelset_geometry;
    void init_eb(const Geometry &geom,const BoxArray &ba,
            const DistributionMapping &dm);
    void make_wedge_hopper_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real get_levelset_value(amrex::Array4<amrex::Real> phi,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
        amrex::Real ptxyz[AMREX_SPACEDIM],
        amrex::Real lsref)
{
    //assuming a redistribute was done and no
    //particle is outside of the tilebox
    //
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> ls_dx;
    const int lo[AMREX_SPACEDIM]={0};
    const int hi[AMREX_SPACEDIM]={0};
    const GpuArray<int, AMREX_SPACEDIM> order_scheme_pseudo={1,1,1};
    const GpuArray<int, AMREX_SPACEDIM> periodic_pseudo={1,1,1};

    ls_dx[XDIR]=dx[XDIR]/lsref;
    ls_dx[YDIR]=dx[YDIR]/lsref;
    ls_dx[ZDIR]=dx[ZDIR]/lsref;

    int i=amrex::Math::floor((ptxyz[XDIR]-problo[XDIR]+TINYVAL)/ls_dx[XDIR]);
    int j=amrex::Math::floor((ptxyz[YDIR]-problo[YDIR]+TINYVAL)/ls_dx[YDIR]);
    int k=amrex::Math::floor((ptxyz[ZDIR]-problo[ZDIR]+TINYVAL)/ls_dx[ZDIR]);

    amrex::Real lsval=0.0;
    for(int n=0;n<2;n++)
    {
        for(int m=0;m<2;m++)
        {
            for(int l=0;l<2;l++)
            {
                lsval += phi(i+l,j+m,k+n)*
                    basisval(l,m,n,i,j,k,ptxyz,problo,ls_dx,order_scheme_pseudo,periodic_pseudo,lo,hi);
            }
        }
    }

    return(lsval);
}

AMREX_GPU_DEVICE 
AMREX_FORCE_INLINE
void get_levelset_grad(amrex::Array4<amrex::Real> phi,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
        amrex::Real ptxyz[AMREX_SPACEDIM],
        amrex::Real lsref,amrex::Real lsgrad[AMREX_SPACEDIM])
{
    //assuming a redistribute was done and no
    //particle is outside of the tilebox
    //
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> ls_dx;
    const int lo[AMREX_SPACEDIM]={0};
    const int hi[AMREX_SPACEDIM]={0};
    
    const GpuArray<int, AMREX_SPACEDIM> order_scheme_pseudo={1,1,1};
    const GpuArray<int, AMREX_SPACEDIM> periodic_pseudo={1,1,1};

    
=======
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Translation.H>

#include <AMReX_EB_utils.H>
#include <constants.H>
#include <interpolants.H>

using namespace amrex;

namespace mpm_ebtools
{
    extern EBFArrayBoxFactory* ebfactory;
    extern MultiFab *lsphi;
    extern int ls_refinement;
    extern const int max_blocks;
    extern bool using_levelset_geometry;
    void init_eb(const Geometry &geom,const BoxArray &ba,
            const DistributionMapping &dm);
    void make_wedge_hopper_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real get_levelset_value(amrex::Array4<amrex::Real> phi,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
        amrex::Real ptxyz[AMREX_SPACEDIM],
        amrex::Real lsref)
{
    //assuming a redistribute was done and no
    //particle is outside of the tilebox
    //
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> ls_dx;
    const int lo[AMREX_SPACEDIM]={0};
    const int hi[AMREX_SPACEDIM]={0};
    const Array<int, AMREX_SPACEDIM> order_scheme_pseudo={1,1,1};
    const Array<int, AMREX_SPACEDIM> periodic_pseudo={1,1,1};

    ls_dx[XDIR]=dx[XDIR]/lsref;
    ls_dx[YDIR]=dx[YDIR]/lsref;
    ls_dx[ZDIR]=dx[ZDIR]/lsref;

    int i=amrex::Math::floor((ptxyz[XDIR]-problo[XDIR]+TINYVAL)/ls_dx[XDIR]);
    int j=amrex::Math::floor((ptxyz[YDIR]-problo[YDIR]+TINYVAL)/ls_dx[YDIR]);
    int k=amrex::Math::floor((ptxyz[ZDIR]-problo[ZDIR]+TINYVAL)/ls_dx[ZDIR]);

    amrex::Real lsval=0.0;
    for(int n=0;n<2;n++)
    {
        for(int m=0;m<2;m++)
        {
            for(int l=0;l<2;l++)
            {
                lsval += phi(i+l,j+m,k+n)*
                    basisval(l,m,n,i,j,k,ptxyz,problo,ls_dx,order_scheme_pseudo,periodic_pseudo,lo,hi);
            }
        }
    }

    return(lsval);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void get_levelset_grad(amrex::Array4<amrex::Real> phi,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
        amrex::Real ptxyz[AMREX_SPACEDIM],
        amrex::Real lsref,amrex::Real lsgrad[AMREX_SPACEDIM])
{
    //assuming a redistribute was done and no
    //particle is outside of the tilebox
    //
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> ls_dx;
    const int lo[AMREX_SPACEDIM]={0};
    const int hi[AMREX_SPACEDIM]={0};
    const Array<int, AMREX_SPACEDIM> order_scheme_pseudo={1,1,1};
    const Array<int, AMREX_SPACEDIM> periodic_pseudo={1,1,1};

>>>>>>> refs/heads/main
    ls_dx[XDIR]=dx[XDIR]/lsref;
    ls_dx[YDIR]=dx[YDIR]/lsref;
    ls_dx[ZDIR]=dx[ZDIR]/lsref;

    int i=amrex::Math::floor((ptxyz[XDIR]-problo[XDIR]+TINYVAL)/ls_dx[XDIR]);
    int j=amrex::Math::floor((ptxyz[YDIR]-problo[YDIR]+TINYVAL)/ls_dx[YDIR]);
    int k=amrex::Math::floor((ptxyz[ZDIR]-problo[ZDIR]+TINYVAL)/ls_dx[ZDIR]);

    lsgrad[XDIR]=0.0;
    lsgrad[YDIR]=0.0;
    lsgrad[ZDIR]=0.0;
    for(int n=0;n<2;n++)
    {
        for(int m=0;m<2;m++)
        {
            for(int l=0;l<2;l++)
            {
                lsgrad[XDIR] += phi(i+l,j+m,k+n)*
                    basisvalder(XDIR,l,m,n,i,j,k,ptxyz,problo,ls_dx,order_scheme_pseudo,periodic_pseudo,lo,hi);
                lsgrad[YDIR] += phi(i+l,j+m,k+n)*
                    basisvalder(YDIR,l,m,n,i,j,k,ptxyz,problo,ls_dx,order_scheme_pseudo,periodic_pseudo,lo,hi);
                lsgrad[ZDIR] += phi(i+l,j+m,k+n)*
                    basisvalder(ZDIR,l,m,n,i,j,k,ptxyz,problo,ls_dx,order_scheme_pseudo,periodic_pseudo,lo,hi);
            }
        }
    }
}
#endif
