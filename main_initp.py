#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import Integrator
import Detector
import Params

Detector.LoadBField("bfield/bfield.txt")

print "\nDone loading B-field...\nB-field at center:", Detector.getBField(0,0,0), '\n'

# define the initial momenta (in MeV)
init_p = []
init_p.append([2000, 0, 0])
#init_p.append([4800, 0, 2000])
init_p.append([10000, 0, 0])
init_p.append([10000, 5000, 17000])
init_p.append([30000, 40000, 5000])

colors = ['r', 'g', 'b', 'c']

print 'Initial Momenta (colors r,g,b,c):'
for i in range(len(init_p)):
    print ' -',round(np.linalg.norm(init_p[i])/1000, 1), "GeV"

dt = 0.2
nsteps = 250

trajs = list(range(len(init_p)))
trajs_noMSC = list(range(len(init_p)))

for i in range(len(init_p)):
    x0 = np.array([0,0,0]+init_p[i])
    trajs[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
    trajs_noMSC[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)

plt.figure(1)

time = np.arange(0,dt*nsteps+1e-10, dt)
for i in range(len(init_p)):
    plt.plot(time, trajs[i][2,:]-trajs_noMSC[i][2,:], color=colors[i])

plt.xlabel('time (ns)')
plt.ylabel('Difference (m)')
plt.title('Difference in z-coordinate vs. t (with/without MSC)')

fig = plt.figure(2, figsize=(14.5,5.5))

## xz slice
x = np.arange(-Params.RMAX, Params.RMAX+1e-10, Params.DR)/100
z = np.arange(Params.ZMIN, Params.ZMAX+1e-10, Params.DZ)/100
    
Z,X = np.meshgrid(z,x)

plt.subplot2grid((1,5),(0,0),colspan=3)

# draw mag field

mag = np.append(Params.Bmag[::-1,:,0],Params.Bmag[1:,:,180/Params.DPHI],0)
bmplot = plt.pcolor(Z,X,mag,cmap='afmhot',vmax=5.0)
#bmcb = plt.colorbar(bmplot, orientation='horizontal')

#if np.linalg.norm(Detector.getBField(0,0,0)) != 0:
#    plt.quiver(Z,X,Bzy,Bxy)

sl = Params.solLength
sr = Params.solRad

# draw trajectory
for i in range(len(init_p)):
    plt.plot(trajs[i][2,:],trajs[i][0,:],'-', linewidth=2, color=colors[i])

plt.axis([-15,15,-9,9])
plt.xlabel("z (m)")
plt.ylabel("x (m)")

## xy slice

plt.subplot2grid((1,5),(0,3), colspan=2)

# draw solenoid outline
t = np.linspace(0, 2*np.pi, 100)
plt.plot(sr*np.cos(t),sr*np.sin(t), '-')

# draw trajectory
for i in range(len(init_p)):
    plt.plot(trajs[i][0,:],trajs[i][1,:],'-', linewidth=2, color=colors[i])

plt.axis([-9,9,-9,9])
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
ax.plot(xs=sr*np.cos(t), ys=sr*np.sin(t), zs=sl/2, color='k')
ax.plot(xs=sr*np.cos(t), ys=sr*np.sin(t), zs=-sl/2, color='k')
for i in range(8):
    th = i * 2*np.pi/8
    x = sr*np.cos(th)
    y = sr*np.sin(th)
    ax.plot(xs=[x,x], ys=[y,y], zs=[-sl/2, sl/2], color='k')

ax.set_xlim(-9,9)
ax.set_ylim(-9,9)
ax.set_zlim(-15,15)
plt.show()



