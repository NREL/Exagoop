import numpy as np
import sys
import glob

class AMReXParticleHeader(object):
    
    def __init__(self, header_filename):
        self.real_component_names = []
        self.int_component_names = []
        with open(header_filename, "r") as f:
            self.version_string = f.readline().strip()

            particle_real_type = self.version_string.split('_')[-1]
            particle_real_type = self.version_string.split('_')[-1]
            if particle_real_type == 'double':
                self.real_type = np.float64
            elif particle_real_type == 'single':
                self.real_type = np.float32
            else:
                raise RuntimeError("Did not recognize particle real type.")
            self.int_type = np.int32

            self.dim = int(f.readline().strip())
            self.num_int_base = 2
            self.num_real_base = self.dim
            self.num_real_extra = int(f.readline().strip())
            for i in range(self.num_real_extra):
                self.real_component_names.append(f.readline().strip())
            self.num_int_extra = int(f.readline().strip())
            for i in range(self.num_int_extra):
                self.int_component_names.append(f.readline().strip())
            self.num_int = self.num_int_base + self.num_int_extra
            self.num_real = self.num_real_base + self.num_real_extra
            self.is_checkpoint = bool(int(f.readline().strip()))
            self.num_particles = int(f.readline().strip())
            self.max_next_id = int(f.readline().strip())
            self.finest_level = int(f.readline().strip())
            self.num_levels = self.finest_level + 1

            if not self.is_checkpoint:
                self.num_int_base = 0
                self.num_int_extra = 0
                self.num_int = 0

            self.grids_per_level = np.zeros(self.num_levels, dtype='int64')
            for level_num in range(self.num_levels):
                self.grids_per_level[level_num] = int(f.readline().strip())

            self.grids = [[] for _ in range(self.num_levels)]
            for level_num in range(self.num_levels):
                for grid_num in range(self.grids_per_level[level_num]):
                    entry = [int(val) for val in f.readline().strip().split()]
                    self.grids[level_num].append(tuple(entry))

def read_amrex_binary_particle_file(fn):

    timestamp=float(fn.split("plt")[1])
    base_fn = fn + "/particles"
    header = AMReXParticleHeader(base_fn + "/Header")

    
    idtype = "(%d,)i4" % header.num_int    
    if header.real_type == np.float64:
        fdtype = "(%d,)f8" % header.num_real
    elif header.real_type == np.float32:
        fdtype = "(%d,)f4" % header.num_real
    
    idata = np.empty((header.num_particles, header.num_int ))
    rdata = np.empty((header.num_particles, header.num_real))
    
    ip = 0
    for lvl, level_grids in enumerate(header.grids):
        for (which, count, where) in level_grids:
            if count == 0: continue
            fn = base_fn + "/Level_%d/DATA_%05d" % (lvl, which)

            with open(fn, 'rb') as f:
                f.seek(where)
                if header.is_checkpoint:
                    ints   = np.fromfile(f, dtype = idtype, count=count)
                    idata[ip:ip+count] = ints

                floats = np.fromfile(f, dtype = fdtype, count=count)
                rdata[ip:ip+count] = floats            
            ip += count

    return idata, rdata, timestamp

def write_paraview_file_particles(fname,pts,ncdata):

    outfile=open(fname,'w')

    Npts=pts.shape[0]

    outfile.write("<?xml version=\"1.0\"?>\n")
    outfile.write("<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
    outfile.write("<PolyData>\n")
    outfile.write("<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"%(Npts))

    outfile.write("<PointData>\n")
    n_ncdata=ncdata.shape[0]
    if(n_ncdata > 0):
        for ndataset in range(n_ncdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"Point_data%d\" format=\"ascii\">\n"%(ndataset))
            for i in range(ncdata.shape[1]):
                if(ncdata[ndataset][i]<1e-50):
                    ncdata[ndataset][i]=0.0
                outfile.write("%e "%(ncdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</PointData>\n")

    outfile.write("<Points>\n")
    outfile.write("<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n")
    for i in range(Npts):
        outfile.write("%e\t%e\t%e\t"%(pts[i][0],pts[i][1],pts[i][2]))
    outfile.write("\n</DataArray>\n")
    outfile.write("</Points>\n")

    outfile.write("</Piece>\n")
    outfile.write("</PolyData>\n")
    outfile.write("</VTKFile>\n")

    outfile.close()


if __name__ == "__main__":
    fn_pattern = sys.argv[1]
    varnum=int(sys.argv[2])
    fn_list = glob.glob(fn_pattern)
    fn_list.sort()

    outfile=open("avg_quantities.dat","w")

    for i, fn in enumerate(fn_list):
        print("reading %s"%(fn))
        idata, rdata, timestamp = read_amrex_binary_particle_file(fn)
        ppos = rdata[:,0:3]   # assumes 3D
        ncdata = np.transpose(rdata[:,3+varnum:3+varnum+1])
        print(ncdata.shape)
        write_paraview_file_particles(fn+".vtp",ppos,ncdata)
