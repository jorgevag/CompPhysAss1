import os,sys
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

f=open('Output/flash.txt', 'r').readlines()
N=len(f)
t = np.zeros(N)
potential = np.zeros(N)
for i in range(0,N):
    w=f[i].split()
    t[i]= w[0]
    potential[i] = w[1]
    #try:
        #x=[float(j) for j in w[0]]
        #potential=[float(j) for j in w[1]]
    #except ValueError,e:
        #print "error",e,"on line",i

fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_xlabel("t")
ax1.set_ylabel("f(t)")
#ax1.axis([-0.5,0.5,0, 2])
ax1.plot(t, potential, c='r', label='f(t)')

leg = ax1.legend()
plt.show()

