import numpy as np
import matplotlib.pyplot as plt
from sys import argv

data=np.loadtxt('AxialBarEnergy.out.0')
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.grid('on')

plt.plot(data[:,0],data[:,1],label='TKE',color='r')
plt.plot(data[:,0],data[:,2],label='TSE',color='blue')
plt.plot(data[:,0],data[:,3],label='TE',color='black')
plt.xlabel("Time ")
plt.ylabel("Energy ")
lgd = ax.legend()  
# saving the file.Make sure you 
# use savefig() before show().
plt.savefig(argv[1])

