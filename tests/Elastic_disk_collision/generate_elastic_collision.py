import numpy as np
from sys import argv
from _locale import ABMON_10

# Use the same parameters as in input file
blo    = np.array([float(0.0),float(0.0),float(0.0)])
bhi    = np.array([float(1.4),float(1.4),float(1.4)])
ncells = np.array([14,14,14])
npart  = 0 
dx = (bhi-blo)/ncells;
if(dx[0]!=dx[1] or dx[0]!=dx[2] or dx[1]!=dx[2]):
    print("Error! mesh sizes are not same in all directions",dx[0],dx[1],dx[2])
nparticle_per_cells_eachdir=1
xmin=0.0
xmax=1.4
ymin=0.0
ymax=1.4
zmin=0
zmax=1.4

print(range(nparticle_per_cells_eachdir))

print('Number of particles = ',npart)
outfile=open("mpm_particles.dat","w")
outfile.write("%d\n"%(npart));

dens=997.5
phase=0
rad=0.025
E = 1000
nu=0.3
#Volume in each cell
xc=0.4
yc=0.4
zc=0.7
vol_cell=dx[0]*dx[1]*dx[2]
vol_particle=vol_cell/(nparticle_per_cells_eachdir*nparticle_per_cells_eachdir*nparticle_per_cells_eachdir)
rad=(3.0/4.0*vol_particle/3.1416)**(1.0/3.0)
npart=0
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
                            
                            if(((cell_cx-xc)*(cell_cx-xc)+(cell_cy-yc)*(cell_cy-yc)+(cell_cz-zc)*(cell_cz-zc))**0.5<=0.2):
                                npart=npart+1
                                velx=0.1;
                                vely=0.1;
                                velz=0.0;
                                
                                outfile.write("%d\t%e\t%e\t%e\t"%(phase,cell_cx,cell_cy,cell_cz));
                                outfile.write("%e\t%e\t"%(rad,dens));
                                outfile.write("%e\t%e\t%e\t"%(velx,vely,velz));                                
                                outfile.write("%d\t%e\t%e\n"%(0,E,nu));
xc=1.0
yc=1.0
zc=0.7
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
                            
                            if(((cell_cx-xc)*(cell_cx-xc)+(cell_cy-yc)*(cell_cy-yc)+(cell_cz-zc)*(cell_cz-zc))**0.5<=0.2):
                                npart=npart+1
                                velx=-0.1;
                                vely=-0.1;
                                velz=0.0;
                                
                                outfile.write("%d\t%e\t%e\t%e\t"%(phase,cell_cx,cell_cy,cell_cz));
                                outfile.write("%e\t%e\t"%(rad,dens));
                                outfile.write("%e\t%e\t%e\t"%(velx,vely,velz));                                
                                outfile.write("%d\t%e\t%e\n"%(0,E,nu));

print(npart)
outfile.close()
