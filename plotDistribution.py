#!/usr/bin/env python
#import os,sys
#from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


# Reading Data
f=open('Output/distribution.txt', 'r')
lines = f.readlines()
f.close()
N=len(lines)-1
x = np.zeros(N)
for i in range(N):
    w=lines[i].split()
    x[i] = w[0]


fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_title("Distribution of Particles")
ax1.set_xlabel("$Position$")
ax1.set_ylabel("$Number in bin$")
#ax1.axis([-0.5,0.5,-10, 10])
n, bins, patches = ax1.hist(x, 500, normed=0, facecolor='g', alpha=0.75, label='Distribution of particles')

leg = ax1.legend()
plt.show()


