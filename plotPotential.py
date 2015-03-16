import os,sys
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

f=open('Output/potential.txt', 'r').readlines()
N=len(f)-1
x = np.zeros(N)
potential = np.zeros(N)
for i in range(0,N):
    w=f[i].split()
    x[i]= w[0]
    potential[i] = w[1]
    #try:
        #x=[float(j) for j in w[0]]
        #potential=[float(j) for j in w[1]]
    #except ValueError,e:
        #print "error",e,"on line",i

fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_title("Potential")
ax1.set_xlabel("x")
ax1.set_ylabel("U(x)")
#ax1.axis([-0.5,0.5,0, 2])
ax1.plot(x, potential, c='r', label='potential')
#ax1.plot(data, 'ro', c='b', label='the data')

leg = ax1.legend()
plt.show()
