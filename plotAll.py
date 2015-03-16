import os,sys
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


f=open('Output/force.txt', 'r')
forceLines = f.readlines()
f.close()
N=len(forceLines)-1
x = np.zeros(N)
originLine = np.zeros(N)
force = np.zeros(N)
for i in range(0,N):
    w=forceLines[i].split()
    x[i]= w[0]
    force[i] = w[1]

f=open('Output/potential.txt', 'r')
potentialLines = f.readlines()
f.close()
potential = np.zeros(N)
for i in range(0,N):
    w=potentialLines[i].split()
    x[i]= w[0]
    potential[i] = w[1]


fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_title("Potential and Force")
ax1.set_xlabel("x")
ax1.set_ylabel("U(x),F(x)")
ax1.axis([-0.5,0.5,-10, 10])
ax1.plot(x, potential, c='g', label='potential')
ax1.plot(x, force, c='r', label='force')
ax1.plot(x, originLine, '--', c='b')
#ax1.plot(data, 'ro', c='b', label='the data')

leg = ax1.legend()
plt.show()
#ax1.plot(data, 'ro', c='b', label='the data')
