import numpy as np
from sys import argv

dx = 0.014    # as in input file
no_of_nodes_in_z = 5

r_min=0.0
r_max=0.2
theta_min=0.0
theta_max=360.0 #in degrees
zmin=0.0+dx/2
zmax=zmin+(no_of_nodes_in_z-1)*dx
no_of_nodes_in_r = 21
no_of_nodes_in_t = 73

dr=(r_max-r_min)/(no_of_nodes_in_r-1.0)
dt=(theta_max-theta_min)/(no_of_nodes_in_t-1.0)

xc=0.4
yc=0.4
zc= zmin

npart=no_of_nodes_in_z*no_of_nodes_in_r*no_of_nodes_in_t*2
print('Number of particles = ',npart)
outfile=open("mpm_particles.dat","w")
outfile.write("%d\n"%(npart))

dens=1000.0
phase=0

print('Mass of each particle ',np.pi*r_max*r_max*dx*2*dens/(npart/2*dens))
rad=(np.pi*r_max*r_max*dx*2*dens/(npart/2*dens)*3.0/4.0/np.pi)**(1.0/3.0)
print('Rad = ',rad)
K_BM=2e4
Gama_Pressure=7.0
Dyn_visc=0.001
E = 1000
nu=0.3

for z in range(no_of_nodes_in_z):
    for r in range(no_of_nodes_in_r):
        rval=r_min+r*dr
        for t in range(no_of_nodes_in_t):
            tval=theta_min+t*dt            
            x=xc+rval*np.cos(tval*np.pi/180)
            y=yc+rval*np.sin(tval*np.pi/180)
            zval=zc+z*dx
            
            velx=0.1;
            vely=0.1;
            velz=0.0;
            
            outfile.write("%d\t%e\t%e\t%e\t"%(phase,x,y,zval))
            outfile.write("%e\t%e\t"%(rad,dens))
            outfile.write("%e\t%e\t%e\t"%(velx,vely,velz))            
            outfile.write("%d\t%e\t%e\n"%(0,E,nu))

xc=1.0
yc=1.0           
for z in range(no_of_nodes_in_z):
    for r in range(no_of_nodes_in_r):
        rval=r_min+r*dr
        for t in range(no_of_nodes_in_t):
            tval=theta_min+t*dt            
            x=xc+rval*np.cos(tval*np.pi/180)
            y=yc+rval*np.sin(tval*np.pi/180)
            zval=zc+z*dx

            velx=-0.1;
            vely=-0.1;
            velz=0.0;

            outfile.write("%d\t%e\t%e\t%e\t"%(phase,x,y,zval))
            outfile.write("%e\t%e\t"%(rad,dens))
            outfile.write("%e\t%e\t%e\t"%(velx,vely,velz))            
            outfile.write("%d\t%e\t%e\n"%(0,E,nu))
    
outfile.close()
