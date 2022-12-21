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
<<<<<<< HEAD
    
    void make_wedge_hopper_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int ls_ref = ls_refinement;
        // Define nGrow of level-set and EB
        int nghost = 1;
        
        const auto plo = geom.ProbLoArray();
        const auto phi = geom.ProbHiArray();

        amrex::Real exit_size     = 0.0002;
        amrex::Real bin_size      = 0.0002;
        amrex::Real funnel_height  = 0.0002;
        amrex::Real vertoffset=0.5*(plo[1]+phi[1]);
        Vector<amrex::Real> centervec(3);

        amrex::ParmParse pp("wedge_hopper");
        pp.get("exit_size",     exit_size);
        pp.get("bin_size",      bin_size);
        pp.get("funnel_height",  funnel_height);
        pp.get("vertical_offset", vertoffset);
        
        Array<amrex::Real,3> funnel_point1 ={0.5*exit_size,0.0,0.0};
        Array<amrex::Real,3> funnel_normal1={funnel_height,0.5*(exit_size-bin_size),0.0};
        EB2::PlaneIF funnel1(funnel_point1, funnel_normal1);
        Array<amrex::Real,3> bin_point1={0.5*bin_size,funnel_height,0.0};
        Array<amrex::Real,3> bin_normal1={1.0,0.0,0.0};
        EB2::PlaneIF bin1(bin_point1, bin_normal1);
        
        Array<amrex::Real,3> funnel_point2 ={-0.5*exit_size,0.0,0.0};
        Array<amrex::Real,3> funnel_normal2={-funnel_height,0.5*(exit_size-bin_size),0.0};
        EB2::PlaneIF funnel2(funnel_point2, funnel_normal2);
        Array<amrex::Real,3> bin_point2={-0.5*bin_size,funnel_height,0.0};
        Array<amrex::Real,3> bin_normal2={-1.0,0.0,0.0};
        EB2::PlaneIF bin2(bin_point2, bin_normal2);

        Array<Real,3> center={0.5*(plo[0]+phi[0]),vertoffset,0.5*(plo[2]+phi[2])};
        auto hopper_alone = EB2::translate(EB2::makeUnion(funnel1,bin1,funnel2,bin2),center);
        
        amrex::Real len[AMREX_SPACEDIM]={phi[0]-plo[0],phi[1]-plo[1],phi[2]-plo[2]};

        RealArray lo,hi;
        lo[0]=plo[0]-len[0];
        lo[1]=plo[1]-len[1];
        lo[2]=plo[2]-len[2];

        hi[0]=phi[0]+len[0];
        hi[1]=vertoffset;
        hi[2]=phi[2]+len[2];
        EB2::BoxIF box_below(lo, hi, false); 

        auto hopper=EB2::makeComplement(EB2::makeUnion(EB2::makeComplement(hopper_alone),box_below));
        
        //Define EB
        //auto hopper_gshop = EB2::makeShop(box_below);
        auto hopper_gshop = EB2::makeShop(hopper);

        //make domain finer for levelset
        Box dom_ls = geom.Domain();
        dom_ls.refine(ls_ref);
        Geometry geom_ls(dom_ls);

        int required_coarsening_level = 0;
        int max_coarsening_level=10;
        if (ls_refinement > 1) 
        {
            int tmp = ls_refinement;
            while (tmp >>= 1) ++required_coarsening_level;
        }

        // Build EB
        EB2::Build(hopper_gshop, geom_ls, 
                   required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
        const EB2::Level &      eblev  = ebis.getLevel(geom);
        //create lslev
        const EB2::Level & lslev = ebis.getLevel(geom_ls);

        //build factory
        ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                                           {nghost, nghost,
                                               nghost}, EBSupport::full);

        //create nodal multifab with level-set refinement
        BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
        ls_ba.refine(ls_ref);
        lsphi = new MultiFab;
        lsphi->define(ls_ba, dm, 1, nghost);

        amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_ref);

    }

    void init_eb(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int nghost = 1;
        std::string geom_type="all_regular";
        amrex::ParmParse pp("eb2");

        pp.query("geom_type",geom_type);
        pp.query("ls_refinement",ls_refinement);

        if(geom_type!="all_regular")
        {
            if(geom_type=="wedge_hopper")
            {
                using_levelset_geometry=true;
                make_wedge_hopper_levelset(geom,ba,dm);
            }
            else
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
=======
    const int max_blocks=16;

    void make_blocks_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int ls_ref = ls_refinement;
        // Define nGrow of level-set and EB
        int nghost = 1;

        const auto plo = geom.ProbLoArray();
        const auto phi = geom.ProbHiArray();

        Array<Array<amrex::Real,3>,max_blocks> lower_corners;
        Array<Array<amrex::Real,3>,max_blocks> upper_corners;

        amrex::Real len_default = std::max(phi[0]-plo[0],
                                           std::max(phi[1]-plo[1],phi[2]-plo[2]));

        Vector<amrex::Real> centervec(3);
        amrex::Real lenfac=20.0;
        centervec[0]=phi[0]+lenfac*(phi[0]-plo[0]);
        centervec[1]=phi[1]+lenfac*(phi[1]-plo[1]);
        centervec[2]=phi[2]+lenfac*(phi[2]-plo[2]);

        for(int i=0;i<max_blocks;i++)
        {
            lower_corners[i]={centervec[0]-len_default, centervec[1]-len_default, centervec[2]-len_default};
            upper_corners[i]={centervec[0]+len_default, centervec[1]+len_default, centervec[2]+len_default};
        }

        int nblocks;
        bool inside=false;
        Vector<amrex::Real> cvec(6);

        amrex::ParmParse pp("blocks");
        pp.query("num_blocks",nblocks);

        for(int i=0;i<nblocks;i++)
        {
            amrex::Print()<<"reading block "<<i<<"\n";

            std::string corners_str=amrex::Concatenate("corners_",i,2); 

            pp.getarr(corners_str.c_str(),    cvec,  0, 6);
            lower_corners[i]={cvec[0], cvec[1], cvec[2]};
            upper_corners[i]={cvec[3], cvec[4], cvec[5]};
        }

        EB2::BoxIF blocks[]={EB2::BoxIF(lower_corners[0], upper_corners[0], inside),
            EB2::BoxIF(lower_corners[1], upper_corners[1], inside),
            EB2::BoxIF(lower_corners[2], upper_corners[2], inside),
            EB2::BoxIF(lower_corners[3], upper_corners[3], inside),
            EB2::BoxIF(lower_corners[4], upper_corners[4], inside),
            EB2::BoxIF(lower_corners[5], upper_corners[5], inside),
            EB2::BoxIF(lower_corners[6], upper_corners[6], inside),
            EB2::BoxIF(lower_corners[7], upper_corners[7], inside),
            EB2::BoxIF(lower_corners[8], upper_corners[8], inside),
            EB2::BoxIF(lower_corners[9], upper_corners[9], inside),
            EB2::BoxIF(lower_corners[10], upper_corners[10], inside),
            EB2::BoxIF(lower_corners[11], upper_corners[11], inside),
            EB2::BoxIF(lower_corners[12], upper_corners[12], inside),
            EB2::BoxIF(lower_corners[13], upper_corners[13], inside),
            EB2::BoxIF(lower_corners[14], upper_corners[14], inside),
            EB2::BoxIF(lower_corners[15], upper_corners[15], inside)};


        auto blocks_union = EB2::makeUnion(blocks[0],blocks[1],blocks[2],blocks[3],
                                           blocks[4],blocks[5],blocks[6],blocks[7],
                                           blocks[8],blocks[9],blocks[10],blocks[11],
                                           blocks[12],blocks[13],blocks[14],blocks[15]);


        //Define EB
        auto blocks_gshop = EB2::makeShop(blocks_union);

        //make domain finer for levelset
        Box dom_ls = geom.Domain();
        dom_ls.refine(ls_ref);
        Geometry geom_ls(dom_ls);

        int required_coarsening_level = 0;
        int max_coarsening_level=10;
        if (ls_refinement > 1) 
        {
            int tmp = ls_refinement;
            while (tmp >>= 1) ++required_coarsening_level;
        }

        // Build EB
        EB2::Build(blocks_gshop, geom_ls, 
                   required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
        const EB2::Level &      eblev  = ebis.getLevel(geom);
        //create lslev
        const EB2::Level & lslev = ebis.getLevel(geom_ls);

        //build factory
        ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                                           {nghost, nghost,
                                               nghost}, EBSupport::full);

        //create nodal multifab with level-set refinement
        BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
        ls_ba.refine(ls_ref);
        lsphi = new MultiFab;
        lsphi->define(ls_ba, dm, 1, nghost);

        amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_ref);

    }

    void make_wedge_hopper_levelset(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int ls_ref = ls_refinement;
        // Define nGrow of level-set and EB
        int nghost = 1;

        const auto plo = geom.ProbLoArray();
        const auto phi = geom.ProbHiArray();

        amrex::Real exit_size     = 0.0002;
        amrex::Real bin_size      = 0.0002;
        amrex::Real funnel_height  = 0.0002;
        amrex::Real vertoffset=0.5*(plo[1]+phi[1]);
        Vector<amrex::Real> centervec(3);

        amrex::ParmParse pp("wedge_hopper");
        pp.get("exit_size",     exit_size);
        pp.get("bin_size",      bin_size);
        pp.get("funnel_height",  funnel_height);
        pp.get("vertical_offset", vertoffset);

        Array<amrex::Real,3> funnel_point1 ={0.5*exit_size,0.0,0.0};
        Array<amrex::Real,3> funnel_normal1={funnel_height,0.5*(exit_size-bin_size),0.0};
        EB2::PlaneIF funnel1(funnel_point1, funnel_normal1);
        Array<amrex::Real,3> bin_point1={0.5*bin_size,funnel_height,0.0};
        Array<amrex::Real,3> bin_normal1={1.0,0.0,0.0};
        EB2::PlaneIF bin1(bin_point1, bin_normal1);

        Array<amrex::Real,3> funnel_point2 ={-0.5*exit_size,0.0,0.0};
        Array<amrex::Real,3> funnel_normal2={-funnel_height,0.5*(exit_size-bin_size),0.0};
        EB2::PlaneIF funnel2(funnel_point2, funnel_normal2);
        Array<amrex::Real,3> bin_point2={-0.5*bin_size,funnel_height,0.0};
        Array<amrex::Real,3> bin_normal2={-1.0,0.0,0.0};
        EB2::PlaneIF bin2(bin_point2, bin_normal2);

        Array<Real,3> center={0.5*(plo[0]+phi[0]),vertoffset,0.5*(plo[2]+phi[2])};
        auto hopper_alone = EB2::translate(EB2::makeUnion(funnel1,bin1,funnel2,bin2),center);

        amrex::Real len[AMREX_SPACEDIM]={phi[0]-plo[0],phi[1]-plo[1],phi[2]-plo[2]};

        RealArray lo,hi;
        lo[0]=plo[0]-len[0];
        lo[1]=plo[1]-len[1];
        lo[2]=plo[2]-len[2];

        hi[0]=phi[0]+len[0];
        hi[1]=vertoffset;
        hi[2]=phi[2]+len[2];
        EB2::BoxIF box_below(lo, hi, false); 

        auto hopper=EB2::makeComplement(EB2::makeUnion(EB2::makeComplement(hopper_alone),box_below));

        //Define EB
        //auto hopper_gshop = EB2::makeShop(box_below);
        auto hopper_gshop = EB2::makeShop(hopper);

        //make domain finer for levelset
        Box dom_ls = geom.Domain();
        dom_ls.refine(ls_ref);
        Geometry geom_ls(dom_ls);

        int required_coarsening_level = 0;
        int max_coarsening_level=10;
        if (ls_refinement > 1) 
        {
            int tmp = ls_refinement;
            while (tmp >>= 1) ++required_coarsening_level;
        }

        // Build EB
        EB2::Build(hopper_gshop, geom_ls, 
                   required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace & ebis   = EB2::IndexSpace::top();
        const EB2::Level &      eblev  = ebis.getLevel(geom);
        //create lslev
        const EB2::Level & lslev = ebis.getLevel(geom_ls);

        //build factory
        ebfactory = new EBFArrayBoxFactory(eblev, geom, ba, dm,
                                           {nghost, nghost,
                                               nghost}, EBSupport::full);

        //create nodal multifab with level-set refinement
        BoxArray ls_ba = amrex::convert(ba, IntVect::TheNodeVector());
        ls_ba.refine(ls_ref);
        lsphi = new MultiFab;
        lsphi->define(ls_ba, dm, 1, nghost);

        amrex::FillSignedDistance (*lsphi,lslev,*ebfactory,ls_ref);

    }

    void init_eb(const Geometry &geom,const BoxArray &ba,const DistributionMapping &dm)
    {
        int nghost = 1;
        std::string geom_type="all_regular";
        amrex::ParmParse pp("mpm");

        pp.query("kind_of_geometry",geom_type);
        pp.query("ls_refinement",ls_refinement);

        if(geom_type!="all_regular")
        {
            if(geom_type=="wedge_hopper")
            {
                using_levelset_geometry=true;
                make_wedge_hopper_levelset(geom,ba,dm);
            }
            else if(geom_type=="blocks")
            {
                using_levelset_geometry=true;
                make_blocks_levelset(geom,ba,dm);
            }
            else if(geom_type=="eb2")
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
            else
            {
                amrex::Abort("embedded geometry not implemented yet\n");   
>>>>>>> refs/heads/main
            }
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
