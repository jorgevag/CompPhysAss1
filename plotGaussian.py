#!/usr/bin/env python
#import os,sys
#from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


# Reading Data
f=open('Output/gaussianRandom.txt', 'r')
lines = f.readlines()
f.close()
N=len(lines)
data = np.zeros(N)
xtheory = np.zeros(N)
ytheory = np.zeros(N)
for i in range(N):
    w=lines[i].split()
    data[i]= w[0]
    xtheory[i] = w[1]
    ytheory[i] = w[2]

# Fit

fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_title("Gaussian Random Number Histogram (Test)$")
ax1.set_xlabel("$X$")
ax1.set_ylabel("$Number in bin$")
#ax1.axis([-0.5,0.5,-10, 10])
ax1.plot(xtheory, ytheory, label = "theory")
n, bins, patches = ax1.hist(data, 500, normed=1, facecolor='g', alpha=0.75, label='gaussian distribution')

leg = ax1.legend()
plt.show()


