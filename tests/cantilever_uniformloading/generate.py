import numpy as np
from sys import argv

prob_lo=np.zeros(3)
prob_hi=np.zeros(3)
ncells=np.array([0,0,0])

infile=open(argv[1],'r')

xmin=float(argv[2])
xmax=float(argv[3])
ymin=float(argv[4])
ymax=float(argv[5])
zmin=float(argv[6])
zmax=float(argv[7])
refine=int(argv[8])

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

ncells=ncells*refine
dx=(prob_hi-prob_lo)/ncells;
npart=0;

outfile=open("mpm_particles.dat","w")
phase=0
dens=1000.0
radmpm=min(dx[0],min(dx[1],dx[2]))*0.5
str=""
for k in range(ncells[2]):
    for j in range(ncells[1]):
        for i in range(ncells[0]):

            x=prob_lo[0]+i*dx[0];
            y=prob_lo[1]+j*dx[1];
            z=prob_lo[2]+k*dx[2];

            if(y>=ymin and y<ymax and z>=zmin and z<zmax and x>=xmin and x<xmax):

                for n3 in range(2):
                    for n2 in range(2):
                        for n1 in range(2):
                            npart+=1
                            partx=x+0.25*dx[0]+n1*0.5*dx[0];
                            party=y+0.25*dx[1]+n2*0.5*dx[1];
                            partz=z+0.25*dx[2]+n3*0.5*dx[2];

                            #outfile.write("%d\t%e\t%e\t%e\t"%(phase,partx,party,partz));
                            #outfile.write("%e\t%e\t"%(radmpm,dens));
                            #outfile.write("%e\t%e\t%e\n"%(0.0,0.0,0.0))
                            str+="%d\t%e\t%e\t%e\t"%(phase,partx,party,partz)
                            str+="%e\t%e\t"%(radmpm,dens)
                            str+="%e\t%e\t%e\n"%(0.0,0.0,0.0)



#outfile.seek(0)
outfile.write("%d\n"%(npart))
outfile.write(str);
outfile.close()
print("npart:",npart)
