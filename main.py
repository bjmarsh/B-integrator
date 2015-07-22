#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import Integrator
import Detector
from Params import *

print "\nB-field at center:", Detector.getBField([0,0,0]), '\n'

# p0 = [15000., 0, 10000.] ## in MeV
# print "Initial p:",np.linalg.norm(p0),"MeV"

# # initial position and momentum
# x0 = np.array([0,0,0]+p0)

### Plot z-crossing as a function of dt

# dtVals = np.arange(0.01, 5.01+1e-10, 0.01)
# zVals = []

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

### plot y as a function of t

# bestTraj = Integrator.rk4(x0, Integrator.traverseBField, 0.1, 1000)
# time = np.arange(0,100+1e-10,.1)
# plt.figure(1)
# plt.plot(time,bestTraj[1,:])

### histogram of z-crossing vals

# zVals = []
# dt = 0.1
# for i in range(200):
#     bestTraj = Integrator.rk4(x0, Integrator.traverseBField, dt, 800)
#     crossInd = -1
#     for i in range(bestTraj.shape[1]-1):
#         if bestTraj[0,i]<20 and bestTraj[0,i+1]>20:
#             crossInd = i
#             break
#     if crossInd==-1:
#         print "Warning: not enough steps at dt =", dt
#         break
#     crossZ = bestTraj[2,crossInd] + (20-bestTraj[0,crossInd])/(bestTraj[0,crossInd+1]-bestTraj[0,crossInd]) * \
#                                     (bestTraj[2,crossInd+1]-bestTraj[2,crossInd])
#     zVals.append(crossZ)
    
# plt.figure(1)
# plt.hist(zVals, bins=20, histtype='stepfilled')

# define the initial momenta (in MeV)
init_p = []
init_p.append([4800, 0, 2000])
init_p.append([10000, 0, 0])
init_p.append([10000, 5000, 17000])
init_p.append([30000, 40000, 5000])

colors = ['r', 'g', 'b', 'c']

print 'Initial Momenta (colors r,g,b,c):'
for i in range(len(init_p)):
    print ' -',round(np.linalg.norm(init_p[i])/1000, 1), "GeV"

dt = 0.1
nsteps = 1000

trajs = list(range(len(init_p)))
trajs_noMSC = list(range(len(init_p)))

for i in range(len(init_p)):
    x0 = np.array([0,0,0]+init_p[i])
    trajs[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, doMSC=True)
    trajs_noMSC[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, doMSC=True)

plt.figure(1)

time = np.arange(0,dt*nsteps+1e-10, dt)
for i in range(len(init_p)):
    plt.plot(time, trajs[i][2,:]-trajs_noMSC[i][2,:], color=colors[i])

plt.xlabel('time (ns)')
plt.ylabel('Difference (m)')
plt.title('Difference in z-coordinate vs. t (with/without MSC)')

fig = plt.figure(2, figsize=(11,5.5))

## xz slice
x = np.linspace(-20,20,21)
z = np.linspace(-30,30,21)

(X,Z) = np.meshgrid(x,z)

Bzy = np.zeros(X.shape)
Bxy = np.zeros(X.shape)

for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        B = Detector.getBField(np.array([X[i,j],0,Z[i,j]]))
        if np.linalg.norm(B)==0:
            continue
        Bxy[i,j] = B[0]/np.linalg.norm(B)
        Bzy[i,j] = B[2]/np.linalg.norm(B)


plt.subplot(1,2,1)
# draw mag field
if BMag != 0:
    plt.quiver(Z,X,Bzy,Bxy)

# draw solenoid outline
plt.plot([solLength/2, solLength/2, -solLength/2, -solLength/2, solLength/2], 
         [-solRad, solRad, solRad, -solRad, -solRad], '-')

# draw trajectory
for i in range(len(init_p)):
    plt.plot(trajs[i][2,:],trajs[i][0,:],'-', linewidth=2, color=colors[i])

plt.axis([-30,30,-20,20])
plt.xlabel("z (m)")
plt.ylabel("x (m)")

## xy slice

plt.subplot(1,2,2)

# draw solenoid outline
t = np.linspace(0, 2*np.pi, 100)
plt.plot(solRad*np.cos(t),solRad*np.sin(t), '-')

# draw trajectory
for i in range(len(init_p)):
    plt.plot(trajs[i][0,:],trajs[i][1,:],'-', linewidth=2, color=colors[i])

plt.axis([-20,20,-20,20])
plt.xlabel("x (m)")
plt.ylabel("y (m)")

## 3d view

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(num=3, figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

for i in range(len(init_p)):
    ax.plot3D(xs=trajs[i][0,:], ys=trajs[i][1,:], zs=trajs[i][2,:], color=colors[i])

ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')

t = np.linspace(0, 2*np.pi, 100)
ax.plot(xs=solRad*np.cos(t), ys=solRad*np.sin(t), zs=solLength/2, color='k')
ax.plot(xs=solRad*np.cos(t), ys=solRad*np.sin(t), zs=-solLength/2, color='k')
for i in range(8):
    th = i * 2*np.pi/8
    x = solRad*np.cos(th)
    y = solRad*np.sin(th)
    ax.plot(xs=[x,x], ys=[y,y], zs=[-solLength/2, solLength/2], color='k')

ax.set_xlim(-20,20)
ax.set_ylim(-20,20)
ax.set_zlim(-30,30)
plt.show()



