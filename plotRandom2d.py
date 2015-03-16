#import os,sys
#from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


f=open('Output/random2d.txt', 'r')
lines = f.readlines()
f.close()
N=len(lines)-1
x = np.zeros(N)
y = np.zeros(N)
for i in range(N):
    w=lines[i].split()
    x[i]= w[0]
    y[i] = w[1]


fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_title("Random Number Test \n $(x,y) = (r_n, r_{n+1})$")
ax1.set_xlabel("$X$")
ax1.set_ylabel("$Y$")
#ax1.axis([-0.5,0.5,-10, 10])
ax1.scatter(x, y, label = "rand")

leg = ax1.legend()
plt.show()

