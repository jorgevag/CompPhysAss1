#!/usr/bin/env python
#import os,sys
#from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


# Reading Data
f=open('Output/trajectory.txt', 'r')
lines = f.readlines()
f.close()
N=len(lines)-1
x = np.zeros(N)
y = np.zeros(N)
for i in range(N):
    w=lines[i].split()
    x[i] = w[0]
    y[i] = w[1]

# Fit

fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_title("Trajectory")
ax1.set_xlabel("Position $x$ [m]")
ax1.set_ylabel("Time $t$ [s]")
#ax1.axis([-0.5,0.5,-10, 10])
ax1.plot(x, y, label = "trajectory")

leg = ax1.legend()
plt.show()



