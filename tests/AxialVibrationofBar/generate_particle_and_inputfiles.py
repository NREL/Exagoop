import numpy as np
from sys import argv
import os

# Use the same parameters as in input file

no_of_cell_in_x=int(argv[1])
L=25.0
buffer =4
dx1=L/no_of_cell_in_x


blo    = np.array([float(0.0),float(-(buffer+0.5)*dx1),float(-(buffer+0.5)*dx1)])
bhi    = np.array([float(L+buffer*dx1),float((buffer+0.5)*dx1),float((buffer+0.5)*dx1)])
ncells = np.array([no_of_cell_in_x+buffer,2*buffer+1,2*buffer+1])
npart  = 0 
dx = (bhi-blo)/ncells;
if(dx[0]!=dx[1] or dx[0]!=dx[2] or dx[1]!=dx[2]):
    print("Error! mesh sizes are not same in all directions")
nparticle_per_cells_eachdirx=int(argv[2])
nparticle_per_cells_eachdiry=1
nparticle_per_cells_eachdirz=1
xmin=0.0
xmax=L
print(dx1)
ymin=-dx1/2.0
ymax=dx1/2.0
zmin=-dx1/2.0
zmax=dx1/2.0

print("Domain extent: ",xmin,ymin,zmin,xmax,ymax,zmax)
print(range(nparticle_per_cells_eachdirx))
print('No of cells = ',ncells)
for k in range(ncells[2]):
    for j in range(ncells[1]):
        for i in range(ncells[0]):
            c_cx=blo[0]+(i)*dx[0]
            c_cy=blo[1]+(j)*dx[1]
            c_cz=blo[2]+(k)*dx[2]
            if(c_cx>=xmin and c_cx<xmax and c_cy>=ymin and c_cy<ymax and c_cz>=zmin and c_cz<zmax):
                for ii in range(nparticle_per_cells_eachdirx):
                    for jj in range(nparticle_per_cells_eachdiry):
                        for kk in range(nparticle_per_cells_eachdirz):
                            npart=npart+1

print('Number of particles = ',npart)
outfile=open("mpm_particles.dat","w")
outfile.write("%d\n"%(npart));

print("Blo = ",blo)
print("Bhi = ",bhi)


dens=1
phase=0
rad=0.025
E = 100
nu=0.0

#Velocity config
n = 1
beta_n = (2*n-1)/2.0*np.pi/L
print("Beta = ",beta_n)
v0=0.1



#Volume in each cell
vol_cell=dx[0]*dx[1]*dx[2]
vol_particle=vol_cell/(nparticle_per_cells_eachdirx*nparticle_per_cells_eachdiry*nparticle_per_cells_eachdirz)
rad=(3.0/4.0*vol_particle/3.1416)**(1.0/3.0)
for k in range(ncells[2]):
    for j in range(ncells[1]):
        for i in range(ncells[0]):
            c_cx=blo[0]+(i)*dx[0]
            c_cy=blo[1]+(j)*dx[1]
            c_cz=blo[2]+(k)*dx[2]
            if(c_cx>=xmin and c_cx<xmax and c_cy>=ymin and c_cy<ymax and c_cz>=zmin and c_cz<zmax):
                for ii in range(nparticle_per_cells_eachdirx):
                    for jj in range(nparticle_per_cells_eachdiry):
                        for kk in range(nparticle_per_cells_eachdirz):
                            
                            cell_cx=c_cx+(2*ii+1)*dx[0]/(2.0*nparticle_per_cells_eachdirx)
                            cell_cy=c_cy+(2*jj+1)*dx[1]/(2.0*nparticle_per_cells_eachdiry)
                            cell_cy=v0*np.sin(beta_n*cell_cx)
                            cell_cz=c_cz+(2*kk+1)*dx[2]/(2.0*nparticle_per_cells_eachdirz)
                            vely=0.0
                            velx=v0*np.sin(beta_n*cell_cx)
                            velz=0.0
    
                            outfile.write("%d\t%e\t%e\t%e\t"%(phase,cell_cx,cell_cy,0.0));
                            outfile.write("%e\t%e\t"%(rad,dens));
                            outfile.write("%e\t%e\t%e\t"%(velx,vely,velz));
                            outfile.write("%d\t%e\t%e\n"%(0,E,nu));


outfile.close()


xmin=0.0
xmax=L+buffer*dx1
ymin=-dx1*(buffer+0.5)
ymax= dx1*(buffer+0.5)
zmin=-dx1*(buffer+0.5)
zmax= dx1*(buffer+0.5)
outfile=open("inputs","w")
outfile.write("#geometry parameters")
outfile.write("\nmpm.prob_lo = "+str(xmin)+" "+str(ymin)+" "+str(zmin))
outfile.write("\nmpm.prob_hi = "+str(xmax)+" "+str(ymax)+" "+str(zmax))
outfile.write("\nmpm.ncells  = "+str(no_of_cell_in_x+buffer)+" "+str(2*buffer+1)+" "+str(2*buffer+1))
outfile.write("\nmpm.max_grid_size = "+str(no_of_cell_in_x+buffer+1))
outfile.write("\nmpm.is_it_periodic = 0  1  1")

outfile.write("\n#timestepping")
outfile.write("\nmpm.fixed_timestep = 0  #1=> fixed time step, 0=> CFL based adaptive time step")
outfile.write("\nmpm.final_time=50")
outfile.write("\nmpm.max_steps=500000")
outfile.write("\nmpm.write_output_time=0.5")
outfile.write("\nmpm.num_redist = 1")
outfile.write("\nmpm.timestep = 0.00005")
outfile.write("\nmpm.screen_output_time = 0.001")
outfile.write("\nmpm.gravity = 0.0 0.0 0.0")
outfile.write("\nmpm.dtmin=1e-2")
outfile.write("\nmpm.CFL=0.1")
outfile.write("\nmpm.order_scheme="+argv[3])
outfile.write("\nmpm.mass_tolerance = 1e-18")
outfile.write("\nmpm.alpha_pic_flip = "+argv[4])
outfile.write("\nmpm.stress_update_scheme= "+argv[5])
outfile.close()

os.system('./mpm3d.gnu.ex inputs')
os.system('python plot_energy.py '+argv[6])
os.system('python plot_vel.py '+argv[7])

os.system('mv AxialBar.out '+argv[8])
os.system('mv Energy.out '+argv[9])

