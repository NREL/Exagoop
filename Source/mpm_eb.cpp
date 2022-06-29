#include <mpm_eb.H>
#include <AMReX_EB_utils.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

namespace mpm_ebtools
{
    EBFArrayBoxFactory* ebfactory=NULL;
    MultiFab* lsphi=NULL;
    int ls_refinement=1;
    bool using_levelset_geometry=false;

    void init_eb(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int nghost = 1;
        std::string geom_type="all_regular";
        amrex::ParmParse pp("eb2");

        pp.query("geom_type",geom_type);
        pp.query("ls_refinement",ls_refinement);

        if(geom_type!="all_regular")
        {
            using_levelset_geometry=true;
            int required_coarsening_level = 0;
            int max_coarsening_level=10;
            if (ls_refinement > 1) 
            {
                int tmp = ls_refinement;
                while (tmp >>= 1) ++required_coarsening_level;
            }
            Box dom_ls = geom.Domain();
            dom_ls.refine(ls_refinement);
            Geometry geom_ls(dom_ls);
            amrex::EB2::Build(geom_ls, required_coarsening_level, max_coarsening_level);

            const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
            const EB2::Level &      eblev  = ebis.getLevel(geom);
            //create lslev
            const EB2::Level & lslev = ebis.getLevel(geom_ls);

            //build factory
            ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                    {nghost, nghost,nghost}, 
                    EBSupport::full);

            //create nodal multifab with level-set refinement
            BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
            ls_ba.refine(ls_refinement);
            lsphi = new MultiFab;
            lsphi->define(ls_ba, dm, 1, nghost);

            //call signed distance
            amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_refinement);
        }
        
        if(using_levelset_geometry)
        { 
            const std::string& pltfile = "ebplt";
            Box dom_ls = geom.Domain();
            dom_ls.refine(ls_refinement);
            Geometry geom_ls(dom_ls);
            BoxArray plot_ba=ba;
            plot_ba.refine(ls_refinement); 
            MultiFab plotmf(plot_ba,dm,lsphi->nComp(),0);
            amrex::average_node_to_cellcenter(plotmf, 0, *lsphi, 0, lsphi->nComp());
            WriteSingleLevelPlotfile(pltfile, plotmf, {"phi"}, geom_ls, 0.0, 0);
        }
    }
}
