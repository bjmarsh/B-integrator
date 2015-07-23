# Integrator.py
# contains routines to perform the numeric integration

import numpy as np
import Detector
from MultipleScatter import *
import Params

def traverseBField(t, x):
    # x is a 6-element vector (x,y,z,px,py,pz)
    # returns dx/dt
    #
    # if B is in Tesla, dt is in ns, p is in units of MeV/c,  then the basic eq is
    # dp/dt = (89.8755) Qv x B,
    
    
    dxdt = np.zeros(6)

    p = x[3:]
    magp = np.linalg.norm(p)
    E = np.sqrt(magp**2 + Params.m**2)
    v = p/E
    dxdt[:3] = v * 2.9979e-1

    B = Detector.getBField(x[:3])

    dxdt[3:] = (89.8755) * Params.Q * np.cross(v,B)

    return dxdt


# 4th order runge-kutta integrator
def rk4(x0, update_func, dt, nsteps, cutoff=None, cutoffaxis=None):
    # x0 is a vector of initial values e.g. (x0,y0,z0,px0,py0,pz0)
    # update func is as in dx/dt = update_func(x,t)
    # return value is an N by nsteps+1 array, where N is the size of x0
    # each column is x at the next time step

    if cutoff!=None and cutoffaxis==None:
        print "Warning: cutoff axis not specified! Not using cutoff"
        cutoff=None

    x0 = np.array(x0)
    x = np.zeros((x0.size, nsteps+1))
    x[:,0] = x0
    t = 0

    # perform the runge-kutta integration
    for i in range(nsteps):
        k1 = update_func(t, x[:,i])
        k2 = update_func(t+dt/2., x[:,i]+dt*k1/2.)
        k3 = update_func(t+dt/2., x[:,i]+dt*k2/2.)
        k4 = update_func(t+dt, x[:,i]+dt*k3)
        dx_Bfield = dt/6. * (k1 + 2*k2 + 2*k3 + k4)

        # add on the effect of MSC if desired
        dx_MS = np.zeros(x0.size)
        if Params.MSCtype.lower()=='pdg':
            dx_MS = multipleScatterPDG(x[:,i], dt)
        elif Params.MSCtype.lower()=='kuhn':
            dx_MS = multipleScatterKuhn(x[:,i], dt)

        x[:,i+1] = x[:,i] + dx_Bfield + dx_MS

        if cutoff!=None and x[cutoffaxis,i+1]>=cutoff:
            return x[:,:i+2]

        t += dt
        
    if cutoff!=None and x[cutoffaxis,-1]<cutoff:
        print "Warning: cutoff not reached!"

    return x
