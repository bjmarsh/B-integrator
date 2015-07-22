import numpy as np
from GetBField import *
from Params import *

def getNormVectors(v):
    # generate and return a random set of orthogonal vectors in the
    # plane orthogonal to v

    v = np.array(v)

    if np.linalg.norm(v) == 0:
        raise ValueError("Can't generate vector orthogonal to 0 vector!")

    # angular coords of v
    thetav = np.arccos(np.dot(v/np.linalg.norm(v), [0,0,1]))
    phiv = np.arctan2(v[1],v[0])
    
    angle = 2*np.pi*np.random.rand()
    random_unit1 = np.array([np.cos(angle), np.sin(angle), 0]).reshape((3,1))
    random_unit2 = np.array([-np.sin(angle), np.cos(angle), 0]).reshape((3,1))

    # matrix for rotation about y axis by thetav
    rotateY = np.array([[ np.cos(thetav), 0, np.sin(thetav)],
                        [              0, 1,              0],
                        [-np.sin(thetav), 0, np.cos(thetav)]])

    # matrix for rotation about z axis by phiv
    rotateZ = np.array([[np.cos(phiv), -np.sin(phiv), 0],
                        [np.sin(phiv),  np.cos(phiv), 0],
                        [           0,             0, 1]])

    random_unit1 = np.dot(rotateZ, np.dot(rotateY, random_unit1))
    random_unit2 = np.dot(rotateZ, np.dot(rotateY, random_unit2))

    random_unit1 = random_unit1.reshape((3,))
    random_unit2 = random_unit2.reshape((3,))

    return (random_unit1, random_unit2)



def multipleScatter(x, dt):
    # given a velocity and timestep, compute a
    # deflection angle and deviation due to multiple scattering
    # Update the position/velocity and return deflection
    
    global m, Q, X0

    p = x[3:]
    magp = np.linalg.norm(p) # must be in MeV
    E = np.sqrt(magp**2 + m**2)
    v = p/E
    beta = np.linalg.norm(v)

    dx = (beta*2.9979e-7) * dt  # in m
    
    ## in the following we generate a random deflection angle theta
    ## and transvers displacement y for the given momentum. This is taken
    ## from the PDG review chapter on the Passage of Particles through Matter

    # rms of projected theta distribution.
    theta0 = 13.6/(beta*magp) * Q * np.sqrt(dx/X0) * (1 + 0.038*np.log(dx/X0))
    
    # correlation coefficient between theta_plane and y_plane
    rho = 0.87
    
    z1 = np.random.normal()
    z2 = np.random.normal()
    yx = z1*dx*theta0 * np.sqrt((1-rho**2)/3) + z2*rho*dx*theta0/np.sqrt(3)
    thetax = z2*theta0
    
    z1 = np.random.normal()
    z2 = np.random.normal()
    yy = z1*dx*theta0 * np.sqrt((1-rho**2)/3) + z2*rho*dx*theta0/np.sqrt(3)
    thetay = z2*theta0
    
    vx, vy = getNormVectors(p)
    
    # transverse displacement
    disp = yx*vx + yy*vy
    
    # deflection in momentum
    defl = np.linalg.norm(p) * (thetax*vx + thetay*vy)

    return np.append(disp, defl)



def traverseBField(t, x):
    # x is a 6-element vector (x,y,z,px,py,pz)
    # returns dx/dt
    #
    # if B is in Tesla, dt is in ns, p is in units of MeV/c,  then the basic eq is
    # dp/dt = (89.8755) Qv x B,
    # where gamma is the typical relativistic factor
    
    global Q, m
    
    dxdt = np.zeros(6)

    p = x[3:]
    magp = np.linalg.norm(p)
    E = np.sqrt(magp**2 + m**2)
    v = p/E
    dxdt[:3] = v * 2.9979e-1

    B = getBField(x[:3])

    dxdt[3:] = (89.8755) * Q * np.cross(v,B)

    return dxdt


# 4th order runge-kutta integrator
def rk4(x0, update_func, dt, nsteps):
    # x0 is a vector of initial values e.g. (x0,y0,z0,px0,py0,pz0)
    # update func is as in dx/dt = update_func(x,t)
    # return value is an N by nsteps+1 array, where N is the size of x0
    # each column is x at the next time step

    x0 = np.array(x0)
    x = np.zeros((x0.size, nsteps+1))
    x[:,0] = x0
    t = 0

    for i in range(nsteps):
        k1 = update_func(t, x[:,i])
        k2 = update_func(t+dt/2., x[:,i]+dt*k1/2.)
        k3 = update_func(t+dt/2., x[:,i]+dt*k2/2.)
        k4 = update_func(t+dt, x[:,i]+dt*k3)
        dx_Bfield = dt/6. * (k1 + 2*k2 + 2*k3 + k4)

        dx_MS = multipleScatter(x[:,i], dt)

        x[:,i+1] = x[:,i] + dx_Bfield + dx_MS
        t += dt
        

    return x
