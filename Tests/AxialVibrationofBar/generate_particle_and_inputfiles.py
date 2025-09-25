import numpy as np
from sys import argv
import os

# Use the same parameters as in input file

no_of_cell_in_x=int(argv[1])
L=25.0
bufferx =4
buffery = int(argv[2])
bufferz = buffery
dx1=L/no_of_cell_in_x
periodic = int(argv[3])

blo    = np.array([float(0.0),float(-(buffery+0.5)*dx1),float(-(bufferz+0.5)*dx1)])
bhi    = np.array([float(L+bufferx*dx1),float((buffery+0.5)*dx1),float((bufferz+0.5)*dx1)])
ncells = np.array([no_of_cell_in_x+bufferx,2*buffery+1,2*bufferz+1])
npart  = 0 
dx = (bhi-blo)/ncells;
if(dx[0]!=dx[1] or dx[0]!=dx[2] or dx[1]!=dx[2]):
    print("Error! mesh sizes are not same in all directions")
nparticle_per_cells_eachdirx=int(argv[4])
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
                            #cell_cy=v0*np.sin(beta_n*cell_cx)
                            vely=0.0
                            velx=v0*np.sin(beta_n*cell_cx)
                            velz=0.0
    
                            outfile.write("%d\t%e\t%e\t%e\t"%(phase,cell_cx,0.0,0.0));
                            outfile.write("%e\t%e\t"%(rad,dens))
                            outfile.write("%e\t%e\t%e\t"%(velx,vely,velz))
                            outfile.write("%d\t%e\t%e\n"%(0,E,nu))


outfile.close()


xmin=0.0
xmax=L+bufferx*dx1
ymin=-dx1*(buffery+0.5)
ymax= dx1*(buffery+0.5)
zmin=-dx1*(bufferz+0.5)
zmax= dx1*(bufferz+0.5)
outfile=open("inputs_test","w")


outfile.write("#geometry parameters")
outfile.write("\nmpm.prob_lo = "+str(xmin)+" "+str(ymin)+" "+str(zmin)+"\t\t\t#Lower corner of physical domain")
outfile.write("\nmpm.prob_hi = "+str(xmax)+" "+str(ymax)+" "+str(zmax)+"\t\t\t#Upper corner of physical domain")
outfile.write("\nmpm.ncells  = "+str(no_of_cell_in_x+bufferx)+" "+str(2*buffery+1)+" "+str(2*buffery+1))
outfile.write("\nmpm.max_grid_size = "+str(no_of_cell_in_x+bufferx+1))
outfile.write("\nmpm.is_it_periodic = 0  "+str(periodic)+"  "+str(periodic))

outfile.write("\n\n#AMR Parameters")
outfile.write("\n#restart_checkfile = \"\"") 

outfile.write("\n\n#Input files")
outfile.write("\nmpm.use_autogen=0")
outfile.write("\nmpm.mincoords_autogen=0.0 0.0 0.0")
outfile.write("\nmpm.maxcoords_autogen=1.0 1.0 1.0")
outfile.write("\nmpm.vel_autogen=0.0 0.0 0.0") 
outfile.write("\nmpm.constmodel_autogen=0")   
outfile.write("\nmpm.dens_autogen=1.0")      
outfile.write("\nmpm.E_autogen=1e6")        
outfile.write("\nmpm.nu_autogen=0.3")      
outfile.write("\nmpm.bulkmod_autogen=2e6")
outfile.write("\nmpm.Gama_pres_autogen=7")
outfile.write("\nmpm.visc_autogen=0.001")
outfile.write("\nmpm.multi_part_per_cell_autogen=1")
outfile.write("\nmpm.particle_file=\"mpm_particles.dat\"")

outfile.write("\n\n#File output parameters")
outfile.write("\nmpm.prefix_particlefilename=\"./Solution/"+argv[13]+"/plt\"")
outfile.write("\nmpm.prefix_gridfilename=\"./Solution/"+argv[13]+"/nplt\"")             
outfile.write("\nmpm.prefix_densityfilename=\"./Solution/"+argv[13]+"/dens\"")         
outfile.write("\nmpm.prefix_checkpointfilename=\"./Solution/"+argv[13]+"/chk\"")      
outfile.write("\nmpm.num_of_digits_in_filenames=6")

outfile.write("\n\n#Simulation run parameters")
outfile.write("\nmpm.final_time=50.0")                    
outfile.write("\nmpm.max_steps=5000000")                
outfile.write("\nmpm.screen_output_time = 0.001")      
outfile.write("\nmpm.write_output_time=0.5")        
outfile.write("\nmpm.num_redist = 1")                

outfile.write("\n\n#Timestepping parameters")
outfile.write("\nmpm.fixed_timestep = 0")           
outfile.write("\nmpm.timestep = 1.0e-5")           
outfile.write("\nmpm.CFL="+argv[8])                    
outfile.write("\nmpm.dt_min_limit=1e-12")        
outfile.write("\nmpm.dt_max_limit=1e+00")       

outfile.write("\n\n#Numerical schemes")
outfile.write("\nmpm.order_scheme="+argv[5])          
outfile.write("\nmpm.alpha_pic_flip = "+argv[6])  
outfile.write("\nmpm.stress_update_scheme= "+argv[7])
outfile.write("\nmpm.mass_tolerance = 1e-18")

outfile.write("\n\n#Physics parameters")
outfile.write("\nmpm.gravity = 0.0 0.0 0.0")
outfile.write("\nmpm.applied_strainrate_time=0.0")
outfile.write("\nmpm.applied_strainrate=0.0")    
outfile.write("\nmpm.external_loads=0")         
outfile.write("\nmpm.force_slab_lo= 0.0 0.0 0.0")
outfile.write("\nmpm.force_slab_hi= 1.0 1.0 1.0")
outfile.write("\nmpm.extforce = 0.0 0.0 0.0")   

outfile.write("\n\n#Diagnostics and Test")
outfile.write("\nmpm.print_diagnostics= 1")    
outfile.write("\nmpm.is_standard_test= 1")    
outfile.write("\nmpm.test_number= 1")        
outfile.write("\nmpm.axial_bar_E= 100")     
outfile.write("\nmpm.axial_bar_rho= 1")    
outfile.write("\nmpm.axial_bar_L= 25.0")  
outfile.write("\nmpm.axial_bar_modenumber= 1")
outfile.write("\nmpm.axial_bar_v0= 0.1")     

outfile.write("\n\n#Boundary conditions")
outfile.write("\nmpm.bc_lower=1 0 0")
outfile.write("\nmpm.bc_upper=2 0 0")

outfile.close()

os.system('./mpm3d.gnu.MPI.CUDA.ex inputs_test')
os.system('python3 plot_energy.py '+argv[9])
os.system('python3 plot_vel.py '+argv[10])

if(os.path.exists('./PostProc/')):
  if(os.path.exists('./PostProc/'+argv[13])):
    os.system('mv '+argv[9]+' '+'./PostProc/'+argv[13]+'/'+argv[9])
    os.system('mv '+argv[10]+' '+'./PostProc/'+argv[13]+'/'+argv[10])
    os.system('mv AxialBarEnergy.out.0 '+'./PostProc/'+argv[13]+'/'+argv[11])
    os.system('mv AxialBarVel.out.0 '+'./PostProc/'+argv[13]+'/'+argv[12])
  else:
    os.system('mkdir ./PostProc/'+argv[13])    
    os.system('mv '+argv[9]+' '+'./PostProc/'+argv[13]+'/'+argv[9])
    os.system('mv '+argv[10]+' '+'./PostProc/'+argv[13]+'/'+argv[10])
    os.system('mv AxialBarEnergy.out.0 '+'./PostProc/'+argv[13]+'/'+argv[11])
    os.system('mv AxialBarVel.out.0 '+'./PostProc/'+argv[13]+'/'+argv[12])
else:
  os.system('mkdir ./PostProc')
  os.system('mkdir ./PostProc/'+argv[13])
  os.system('mv '+argv[9]+' '+'./PostProc/'+argv[13]+'/'+argv[9])
  os.system('mv '+argv[10]+' '+'./PostProc/'+argv[13]+'/'+argv[10])
  os.system('mv AxialBarEnergy.out.0 '+'./PostProc/'+argv[13]+'/'+argv[11])
  os.system('mv AxialBarVel.out.0 '+'./PostProc/'+argv[13]+'/'+argv[12])

