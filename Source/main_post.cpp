#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <mpm_particle_container.H>
#include <core/Definitions.h>
#include <Partio.h>

using namespace TGSL;

inline void WriteBGEO(const TVP& x, const size_t N, const std::string& filename) {
  Partio::ParticlesDataMutable* parts = Partio::create();
  Partio::ParticleAttribute p_att = parts->addAttribute("position", Partio::VECTOR, 3);
  for (size_t i = 0; i < N; ++i) {
    int idx = parts->addParticle();
    float* p_in = parts->dataWrite<float>(p_att, idx);

    for (int j = 0; j < 3; ++j) {
      p_in[j] = x[i][j];
    }
  }
  Partio::write(filename.c_str(), *parts);
  parts->release();
}

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {

    	//Reading input files for the simulation
        MPMspecs specs;
        specs.read_mpm_specs();	//Read input file

        //
        int coord = 0; //cartesian
        RealBox real_box;
        for (int n = 0; n < AMREX_SPACEDIM; n++)
        {
            real_box.setLo(n, specs.plo[n]);
            real_box.setHi(n, specs.phi[n]);
        }

        IntVect domain_lo(AMREX_D_DECL(0,0,0));
        IntVect domain_hi(AMREX_D_DECL(specs.ncells[XDIR]-1,
                    specs.ncells[YDIR]-1,specs.ncells[ZDIR]-1));
        const Box domain(domain_lo, domain_hi);

        Geometry geom(domain, &real_box, coord, specs.periodic.data());

        BoxArray ba(domain);
        ba.maxSize(specs.max_grid_size);
        DistributionMapping dm(ba);

        int ng_cells = 1;
        if(specs.order_scheme==3)
        {
        	ng_cells = 3;
        }

        //Initialise particle properties
        MPMParticleContainer mpm_pc(geom, dm, ba, ng_cells);
        ParmParse pp;
        std::string restart_file, particle_file="particles";
        pp.get("restart",restart_file);
        bool is_checkpoint = true;
        mpm_pc.Restart(restart_file,particle_file,is_checkpoint);

        const int lev = 0;
        auto& plev  = mpm_pc.GetParticles(lev);
        Long N = 0;

        for(MFIter mfi = mpm_pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
          auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
          N += plev[index].GetArrayOfStructs().numRealParticles();
        }
        ParallelDescriptor::ReduceLongSum(N);
        Print() << "Total np: " << N << std::endl;

        TVP x(N);
        Long cnt = 0;
        for(MFIter mfi = mpm_pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
          const amrex::Box& box = mfi.tilebox();
          int gid = mfi.index();
          int tid = mfi.LocalTileIndex();
          auto index = std::make_pair(gid, tid);
          
          auto& ptile = plev[index];
          auto& aos   = ptile.GetArrayOfStructs();
          const int np = aos.numRealParticles();
          
          auto* pstruct = aos().dataPtr();

          for (int i=0; i<np; ++i)
          {
            const auto& p = pstruct[i];
            x[cnt+i] = {p.pos(XDIR), p.pos(YDIR), p.pos(ZDIR)};
          }
          cnt+=np;
        }

        std::string outfile = restart_file+".bgeo";
        pp.query("outfile",outfile);
        WriteBGEO(x,N,outfile);
    }

    amrex::Finalize();
}
