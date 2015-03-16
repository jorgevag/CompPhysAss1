import os,sys
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

f=open('Output/force.txt', 'r')
lines = f.readlines()
f.close()
N=len(lines)-1
originLine = np.zeros(N)
x = np.zeros(N)
force = np.zeros(N)
for i in range(N):
    w=lines[i].split()
    x[i]= w[0]
    force[i] = w[1]
    #try:
        #x=[float(j) for j in w[0]]
        #potential=[float(j) for j in w[1]]
    #except ValueError,e:
        #print "error",e,"on line",i


fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_title("Force")
ax1.set_xlabel("x")
ax1.set_ylabel("F(x)")
#ax1.axis([-0.5,0.5,-10, 10])
ax1.plot(x, force, c='r', label='force')
ax1.plot(x, originLine, '--', c='b', label='force')
#ax1.plot(data, 'ro', c='b', label='the data')

leg = ax1.legend()
plt.show()

