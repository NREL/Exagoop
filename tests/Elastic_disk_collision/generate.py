import numpy as np
from sys import argv

# Use the same parameters as in input file
blo    = np.array([float(0.0),float(0.0),float(0.0)])
bhi    = np.array([float(1.4),float(1.4),float(1.4)])
ncells = np.array([100,100,100])
npart  = 0 
dx = (bhi-blo)/ncells;
if(dx[0]!=dx[1] or dx[0]!=dx[2] or dx[1]!=dx[2]):
    print("Error! mesh sizes are not same in all directions")
nparticle_per_cells_eachdir=1
xmin=0.0
xmax=1.4
ymin=0.0
ymax=1.4
zmin=0
zmax=1.4

print(range(nparticle_per_cells_eachdir))

for k in range(ncells[2]):
    for j in range(ncells[1]):
        for i in range(ncells[0]):
            c_cx=blo[0]+(i)*dx[0]
            c_cy=blo[1]+(j)*dx[1]
            c_cz=blo[2]+(k)*dx[2]
            if(c_cx>=xmin and c_cx<xmax and c_cy>=ymin and c_cy<ymax and c_cz>=zmin and c_cz<zmax):
                for ii in range(nparticle_per_cells_eachdir):
                    for jj in range(nparticle_per_cells_eachdir):
                        for kk in range(nparticle_per_cells_eachdir):
                            npart=npart+1

print('Number of particles = ',npart)
outfile=open("mpm_particles.dat","w")
outfile.write("%d\n"%(npart));

dens=997.5
phase=0
rad=0.025
K_BM=2e4
Gama_Pressure=7.0
Dyn_visc=0.001
E = 78200000000
nu=0.3
#Volume in each cell
vol_cell=dx[0]*dx[1]*dx[2]
vol_particle=vol_cell/(nparticle_per_cells_eachdir*nparticle_per_cells_eachdir*nparticle_per_cells_eachdir)
rad=(3.0/4.0*vol_particle/3.1416)**(1.0/3.0)
for k in range(ncells[2]):
    for j in range(ncells[1]):
        for i in range(ncells[0]):
            c_cx=blo[0]+(i)*dx[0]
            c_cy=blo[1]+(j)*dx[1]
            c_cz=blo[2]+(k)*dx[2]
            if(c_cx>=xmin and c_cx<xmax and c_cy>=ymin and c_cy<ymax and c_cz>=zmin and c_cz<zmax):
                for ii in range(nparticle_per_cells_eachdir):
                    for jj in range(nparticle_per_cells_eachdir):
                        for kk in range(nparticle_per_cells_eachdir):
                            
                            cell_cx=c_cx+(2*ii+1)*dx[0]/(2.0*nparticle_per_cells_eachdir)
                            cell_cy=c_cy+(2*jj+1)*dx[1]/(2.0*nparticle_per_cells_eachdir)
                            cell_cz=c_cz+(2*kk+1)*dx[2]/(2.0*nparticle_per_cells_eachdir)
                            velx=np.sin(2.0*np.pi/1.4*cell_cx);
                            vely=np.sin(2.0*np.pi/1.4*cell_cy);
                            velz=np.sin(2.0*np.pi/1.4*cell_cz);
    
                            outfile.write("%d\t%e\t%e\t%e\t"%(phase,cell_cx,cell_cy,cell_cz));
                            outfile.write("%e\t%e\t"%(rad,dens));
                            outfile.write("%e\t%e\t%e\t"%(velx,vely,velz));
                            outfile.write("%d\t%e\t%e\n"%(0,E,nu));


outfile.close()
