#!/usr/bin/env python
#import os,sys
#from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


# Reading Data
f=open('Output/energydistribution.txt', 'r')
lines = f.readlines()
f.close()
N=len(lines)
r = np.zeros(N)
x = np.zeros(N)
y = np.zeros(N)
for i in range(N):
    w=lines[i].split()
    r[i] = w[0]
    x[i] = w[1]
    y[i] = w[2]


fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_title("Distribution of Particles")
ax1.set_xlabel("$Position$")
ax1.set_ylabel("$Number in bin$")
#ax1.axis([-0.5,0.5,-10, 10])
n, bins, patches = ax1.hist(r, 500, normed=1, facecolor='g', alpha=0.75, label='Distribution of particles')
ax1.plot(x,y, label='Boltzmann Distribution')

leg = ax1.legend()
plt.show()


