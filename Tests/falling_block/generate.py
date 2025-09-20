import numpy as np
from sys import argv
import sys

def transform(x,y,z,center,angle):
    xd=(x-center[0])*np.cos(angle)+(y-center[1])*np.sin(angle)+center[0]
    yd=-(x-center[0])*np.sin(angle)+(y-center[1])*np.cos(angle)+center[1]
    zd=z
    return(xd,yd,zd)

if(len(argv)<14):
    print("\n\n")
    print("This script generates particles in a box")
    print("****************************************")
    print("Inputs include: box lower corner (3 coordinates), box upper corner (3 coordinates)")
    print("box number of cells (3 ints), inital velocity (3 components), box tilt (in degrees)")
    print("e.g. usage python generate.py 0.4 0.4 0.4 0.6 0.8 0.6 20 20 20 0 -10 0 45")
    print("Hard coded properties: density=500 kg/m3; Youngs modulus=2e5 Pa;")
    print("Poissons ratio:0.2; Particle radius=0.025 m")
    print("\n\n")
    sys.exit()

blo    = np.array([float(argv[1]),float(argv[2]),float(argv[3])])
bhi    = np.array([float(argv[4]),float(argv[5]),float(argv[6])])
ncells = np.array([int(argv[7]),int(argv[8]),int(argv[9])])
angle=float(argv[13])*np.pi/180.0
npart  = ncells[0]*ncells[1]*ncells[2];
center=0.5*(blo+bhi)
dx = (bhi-blo)/ncells;

outfile=open("mpm_particles.dat","w")
outfile.write("%d\n"%(npart));

dens=500
rad=0.025
phase=0
constmodel=0
Ymod=2e5
pratio=0.2
for k in range(ncells[2]):
    for j in range(ncells[1]):
        for i in range(ncells[0]):

            cell_cx=blo[0]+(i+0.5)*dx[0]
            cell_cy=blo[1]+(j+0.5)*dx[1]
            cell_cz=blo[2]+(k+0.5)*dx[2]

            part_cx,part_cy,part_cz=transform(cell_cx,cell_cy,cell_cz,center,angle)

            velx=float(argv[10]);
            vely=float(argv[11]);
            velz=float(argv[12]);
    
            outfile.write("%d\t%e\t%e\t%e\t"%(phase,part_cx,part_cy,part_cz));
            outfile.write("%e\t%e\t"%(rad,dens));
            outfile.write("%e\t%e\t%e\t"%(velx,vely,velz));
            outfile.write("%d\t%e\t%e\n"%(constmodel,Ymod,pratio));


outfile.close()
