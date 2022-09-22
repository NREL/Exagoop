import numpy as np
import matplotlib.pyplot as plt
from sys import argv

data=np.loadtxt('DamBreakWaterfront.out.0')
exp=np.loadtxt('ExperimentalData.dat')

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.grid('on')
id=0
for a in data[:,0]:
  if(a>1.5):
    index = id
    break
  id=id+1
print("index=",id,data[:,0:5])

plt.plot(data[1:id,0],data[1:id,1],label='TKE',color='r')
plt.scatter(exp[:,0],exp[:,1])
plt.xlabel("Time ")
plt.ylabel("X* ")
# saving the file.Make sure you 
# use savefig() before show().
plt.savefig('Waterfront.png')
plt.show()
