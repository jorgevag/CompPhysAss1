#import os,sys
#from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


f=open('Output/random3d.txt', 'r')
lines = f.readlines()
f.close()
N=len(lines)-1
x = np.zeros(N)
y = np.zeros(N)
z = np.zeros(N)
for i in range(N):
    w=lines[i].split()
    x[i]= w[0]
    y[i] = w[1]
    z[i] = w[2]


fig = plt.figure(1)

ax1 = fig.add_subplot(111, projection='3d')
ax1.set_title("Random Number Test \n $(x,y,z) = (r_n, r_{n+1}, r_{n+2})$")
ax1.set_xlabel("$X$")
ax1.set_ylabel("$Y$")
ax1.set_zlabel("$Z$")
#ax1.axis([-0.5,0.5,-10, 10])
#ax1.scatter(x, y)
ax1.scatter(x, y, z, zdir=u'z')

leg = ax1.legend()
plt.show()
