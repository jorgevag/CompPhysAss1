#!/usr/bin/env python
#import os,sys
#from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


# Reading Data
f=open('Output/avdriftvelocity.txt', 'r')
lines = f.readlines()
f.close()
N=len(lines)
t = np.zeros(N)
v = np.zeros(N)
for i in range(N):
    w=lines[i].split()
    t[i] = w[0]
    v[i] = w[1]


fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_xlabel("tau")
ax1.set_ylabel("# of particles/bin")
#ax1.axis([-0.5,0.5,-10, 10])
ax1.plot(t,v,'-o')

#leg = ax1.legend()
plt.show()



