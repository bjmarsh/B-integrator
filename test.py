#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math

def getBField(r):

    global solLength, solRad

    r = np.array(r)
    m = np.array([0.,0.,1500.])
    
    isInsideSol = np.sqrt(r[0]**2+r[1]**2)<=solRad and abs(r[2])<solLength/2

    if abs(r[2])<solLength/2:
        r = np.array([solRad,0,0])
    elif r[2]<-solLength/2:
        r[2] = r[2]+solLength/2
    elif r[2]>solLength/2:
        r[2] = r[2]-solLength/2

    mag = np.linalg.norm(r)
    # if mag<solRad:
    #     mag = solRad
        
    B = 3.0*r*np.dot(m,r)/mag**5 - m/mag**3

    if isInsideSol:
        B = -B

    return B



def update(t, x):
    # x is a 6-element vector (x,y,z,vx,vy,vz)
    # returns dx/dt
    #
    # if B is in Tesla, dt is in ns, v is in units of c, m is in MeV,  then the basic eq is
    # dv/dt = (gamma + v^2*gamma^3)^-1 * (89.875) Qv/m x B,
    # where gamma is the typical relativistic factor
    
    global Q,m

    dxdt = np.zeros(6)
    dxdt[0] = x[3] * 2.998e-1
    dxdt[1] = x[4] * 2.998e-1
    dxdt[2] = x[5] * 2.998e-1

    B = getBField([x[0],x[1],x[2]])

#    print np.linalg.norm(B)  
    if math.isnan(np.linalg.norm(B)):
        exit(0)

    v = np.array([x[3],x[4],x[5]])
    magv = np.linalg.norm(v)
    gamma = 1./(1-magv**2)**0.5
    mult = 1./(gamma+magv**2*gamma**3)

    dxdt[3:] = mult * (89.875) * Q/m * np.cross(v,B)

    return dxdt


def rk4(x0, update_func, dt, nsteps):
    # x0 is a vector of initial values e.g. (x0,y0,z0,vx0,vy0,vz0)
    # update func is as in dx/dt = update_func(x,t)
    # return value is an N by nsteps+1 array, where N is the size of x0
    # each column is x at the next time step

    x0 = np.array(x0)
    x = np.zeros((x0.size, nsteps+1))
    x[:,0] = x0
    t = 0

    for i in range(nsteps):
        k1 = update_func(t, x[:,i])
        k2 = update_func(t+dt/2., x[:,i]+k1/2.)
        k3 = update_func(t+dt/2., x[:,i]+k2/2.)
        k4 = update_func(t+dt, x[:,i]+k3)
        t += dt
        x[:,i+1] = x[:,i]+dt/6.*(k1+2*k2+2*k3+k4)
        print np.linalg.norm(x[3:,i+1])

    return x

solLength = 21.6  ## in m
solRad = 7.3   ## in m

dt = 0.05  ## measured in ns
m = 0.5109989  ## measured in MeV
Q = 1.0  ## measured in units of e

# compute trajectory
x0 = np.array([0,0,0,0.8,0,-0.5998])
traj = rk4(x0, update, dt, 100)

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

# draw mag field
plt.quiver(Z,X,Bzy,Bxy)

# draw solenoid outline
plt.plot([solLength/2, solLength/2, -solLength/2, -solLength/2, solLength/2], [-solRad, solRad, solRad, -solRad, -solRad], '-')

# draw trajectory
plt.plot(traj[2,:],traj[0,:],'-', linewidth=2)

plt.axis([-30,30,-20,20])

## xy slice
x = np.linspace(-20,20,25)
y = np.linspace(-20,20,25)

(X,Y) = np.meshgrid(x,z)

Bxz = np.zeros(X.shape)
Byz = np.zeros(X.shape)

for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        B = getBField(np.array([X[i,j],Y[i,j],0]))
        Bxz[i,j] = B[0]#/np.linalg.norm(B)
        Byz[i,j] = B[1]#/np.linalg.norm(B)

plt.figure(2)

# draw mag field
#plt.quiver(X,Y,Bxz,Byz)

# draw solenoid outline
t = np.linspace(0, 2*np.pi, 100)
plt.plot(solRad*np.cos(t),solRad*np.sin(t), '-')

# draw trajectory
plt.plot(traj[0,:],traj[1,:],'-', linewidth=2)

plt.axis([-20,20,-20,20])


plt.show()



