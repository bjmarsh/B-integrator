#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import Integrator
from GetBField import *
from Params import *


p0 = [20000., 0, 5000] ## in MeV
print "Initial p:",np.linalg.norm(p0),"MeV"

# compute trajectory
x0 = np.array([0,0,0]+p0)
dtVals = np.arange(0.01, 10.01+1e-10, 0.02)
zVals = []

# for dt in dtVals:
#     # print "dt =", dt
#     nsteps = int(10000 * .01/dt)
#     traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
#     #if dt==.01:
#     bestTraj=traj.copy()
#     crossInd = -1
#     for i in range(traj.shape[1]-1):
#         if traj[0,i]<20 and traj[0,i+1]>20:
#             crossInd = i
#             break
#     if crossInd==-1:
#         print "Warning: not enough steps at dt =", dt
#         break
#     crossZ = traj[2,crossInd] + (20-traj[0,crossInd])/(traj[0,crossInd+1]-traj[0,crossInd]) * \
#                                 (traj[2,crossInd+1]-traj[2,crossInd])
#     zVals.append(crossZ)

# plt.figure(1)
# plt.plot(dtVals,zVals)
# plt.xlabel(r"dt")
# plt.ylabel(r"z")

bestTraj = Integrator.rk4(x0, Integrator.traverseBField, 0.1, 1000)
time = np.arange(0,100+1e-10,.1)
plt.figure(1)
plt.plot(time,bestTraj[1,:])

plt.figure(2)

## xz slice
x = np.linspace(-20,20,25)
z = np.linspace(-30,30,25)

(X,Z) = np.meshgrid(x,z)

Bzy = np.zeros(X.shape)
Bxy = np.zeros(X.shape)

for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        B = getBField(np.array([X[i,j],0,Z[i,j]]))
        Bxy[i,j] = B[0]/np.linalg.norm(B)
        Bzy[i,j] = B[2]/np.linalg.norm(B)


plt.subplot(1,2,1)
# draw mag field
plt.quiver(Z,X,Bzy,Bxy)

# draw solenoid outline
plt.plot([solLength/2, solLength/2, -solLength/2, -solLength/2, solLength/2], 
         [-solRad, solRad, solRad, -solRad, -solRad], '-')

# draw trajectory
plt.plot(bestTraj[2,:],bestTraj[0,:],'-', linewidth=2)

plt.axis([-30,30,-20,20])
plt.xlabel("z")
plt.ylabel("x")



## xy slice

plt.subplot(1,2,2)

# draw solenoid outline
t = np.linspace(0, 2*np.pi, 100)
plt.plot(solRad*np.cos(t),solRad*np.sin(t), '-')

# draw trajectory
plt.plot(bestTraj[0,:],bestTraj[1,:],'-', linewidth=2)

plt.axis([-20,20,-20,20])
plt.xlabel("x")
plt.ylabel("y")

plt.show()



