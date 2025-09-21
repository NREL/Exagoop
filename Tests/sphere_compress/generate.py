import numpy as np
from sys import argv

prob_lo=np.zeros(3)
prob_hi=np.zeros(3)
ncells=np.array([0,0,0])

infile=open(argv[1],'r')

ncirc=int(argv[2])
refine=float(argv[3])

for line in infile:
    
    spltline=line.split()

    if(len(spltline)>0):
        if(spltline[0]=="mpm.prob_lo"):
            prob_lo[0]=float(spltline[2])
            prob_lo[1]=float(spltline[3])
            prob_lo[2]=float(spltline[4])

        if(spltline[0]=="mpm.prob_hi"):
            prob_hi[0]=float(spltline[2])
            prob_hi[1]=float(spltline[3])
            prob_hi[2]=float(spltline[4])

        if(spltline[0]=="mpm.ncells"):
            ncells[0]=int(spltline[2])
            ncells[1]=int(spltline[3])
            ncells[2]=int(spltline[4])

infile.close()

ncells=ncells*int(refine)
dx=(prob_hi-prob_lo)/ncells;
mindx=np.min(dx)
prob_len=prob_hi-prob_lo
circradii=0.5*prob_len[0]/ncirc
rad2=circradii**2
npart=0
xcen=np.linspace(prob_lo[0]+circradii,prob_hi[0]-circradii,ncirc)
ycen=np.zeros(ncirc)+circradii
zcen=np.zeros(ncirc)+0.5*(prob_hi[2]+prob_lo[2])


outfile=open("mpm_particles.dat","w")
phase=0
dens=1000.0
str=""
for k in range(ncells[2]):
    for j in range(ncells[1]):
        for i in range(ncells[0]):

            x=prob_lo[0]+(i+0.5)*dx[0];
            y=prob_lo[1]+(j+0.5)*dx[1];
            z=prob_lo[2]+(k+0.5)*dx[2];

            inside=0;
            for c in range(ncirc):
                dist2=(x-xcen[c])**2+(y-ycen[c])**2
                if(dist2<rad2):
                    inside=1
                    break

            if(inside==1):
                npart+=1
                str+="%d\t%e\t%e\t%e\t"%(phase,x,y,z)
                str+="%e\t%e\t"%(mindx/2,dens)
                str+="%e\t%e\t%e\n"%(0.0,0.0,0.0)






#outfile.seek(0)
outfile.write("%d\n"%(npart))
outfile.write(str);
outfile.close()
print("npart:",npart)
