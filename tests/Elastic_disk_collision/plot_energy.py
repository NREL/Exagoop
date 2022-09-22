import numpy as np
import matplotlib.pyplot as plt
data=np.loadtxt('Energy.out')



fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.grid('on')

plt.plot(data[:,0],data[:,1],label='TKE')
plt.plot(data[:,0],data[:,2],label='TSE')
plt.plot(data[:,0],data[:,3],label='TE')
plt.xlabel("Time (s)")
plt.ylabel("Energy (Nm)")
lgd = ax.legend()  
# saving the file.Make sure you 
# use savefig() before show().
plt.savefig("Energy_vs_time_alpha=0.01.png")
plt.show()