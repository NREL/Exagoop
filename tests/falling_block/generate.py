import numpy as np
from sys import argv


blo    = np.array([float(argv[1]),float(argv[2]),float(argv[3])])
bhi    = np.array([float(argv[4]),float(argv[5]),float(argv[6])])
ncells = np.array([int(argv[7]),int(argv[8]),int(argv[9])])
npart  = ncells[0]*ncells[1]*ncells[2];

dx = (bhi-blo)/ncells;
print(dx)

outfile=open("mpm_particles.dat","w")
outfile.write("%d\n"%(npart));

dens=500
rad=0.025
phase=0
for k in range(ncells[2]):
    for j in range(ncells[1]):
        for i in range(ncells[0]):

            cell_cx=blo[0]+(i+0.5)*dx[0]
            cell_cy=blo[1]+(j+0.5)*dx[1]
            cell_cz=blo[2]+(k+0.5)*dx[2]

            velx=0.0;
            vely=-5.0;
            velz=0.0;
    
            outfile.write("%d\t%e\t%e\t%e\t"%(phase,cell_cx,cell_cy,cell_cz));
            outfile.write("%e\t%e\t"%(rad,dens));
            outfile.write("%e\t%e\t%e\n"%(velx,vely,velz));


outfile.close()
