import numpy as np

def update(t, x, Q=1, m=0.511):
    # x is a 6-element vector (x,y,z,px,py,pz)
    # returns dx/dt
    #
    # if B is in Tesla, dt is in ns, p is in units of MeV/c,  then the basic eq is
    # dp/dt = (89.87551787368) Qv x B,
    # where gamma is the typical relativistic factor
    
    dxdt = np.zeros(6)

    p = x[3:]
    magp = np.linalg.norm(p)
    E = np.sqrt(magp**2 + m**2)
    v = p/E
    dxdt[:3] = v * 2.99792458e-1

    B = getBField(x[:3])

#    print np.linalg.norm(B)  
    if math.isnan(np.linalg.norm(B)):
        #exit(0)
        pass

    dxdt[3:] = (89.87551787368) * Q * np.cross(v,B)

    return dxdt



def rk4(x0, update_func, dt, nsteps, Q=1, m=0.511):
    # x0 is a vector of initial values e.g. (x0,y0,z0,px0,py0,pz0)
    # update func is as in dx/dt = update_func(x,t)
    # return value is an N by nsteps+1 array, where N is the size of x0
    # each column is x at the next time step

    x0 = np.array(x0)
    x = np.zeros((x0.size, nsteps+1))
    x[:,0] = x0
    t = 0

    for i in range(nsteps):
        k1 = update_func(t, x[:,i], Q, m)
        k2 = update_func(t+dt/2., x[:,i]+dt*k1/2., Q, m)
        k3 = update_func(t+dt/2., x[:,i]+dt*k2/2., Q, m)
        k4 = update_func(t+dt, x[:,i]+dt*k3, Q, m)
        x[:,i+1] = x[:,i]+dt/6.*(k1+2*k2+2*k3+k4)
        t += dt
        
    return x
