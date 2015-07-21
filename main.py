#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
from GetBField import *
import Integrator



solLength = 21.6  ## in m
solRad = 7.3   ## in m

dt = 0.1  ## measured in ns
m = 0.5109989  ## measured in MeV
Q = 1.0  ## measured in units of e
p0 = [15000., 0, 5000.] ## in MeV
print "Initial p:",np.linalg.norm(p0),"MeV"

# compute trajectory
x0 = np.array([0,0,0]+p0)
traj = Integrator.rk4(x0, update, dt, 1000, Q, m)

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

plt.figure(1)

plt.subplot(1,2,1)
# draw mag field
plt.quiver(Z,X,Bzy,Bxy)

# draw solenoid outline
plt.plot([solLength/2, solLength/2, -solLength/2, -solLength/2, solLength/2], 
         [-solRad, solRad, solRad, -solRad, -solRad], '-')

# draw trajectory
plt.plot(traj[2,:],traj[0,:],'-', linewidth=2)

plt.axis([-30,30,-20,20])




## xy slice

plt.subplot(1,2,2)

# draw solenoid outline
t = np.linspace(0, 2*np.pi, 100)
plt.plot(solRad*np.cos(t),solRad*np.sin(t), '-')

# draw trajectory
plt.plot(traj[0,:],traj[1,:],'-', linewidth=2)

plt.axis([-20,20,-20,20])


plt.show()



