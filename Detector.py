## Detector.py
## methods relating to detector and environment properties

import numpy as np
import Params

def getMaterial(r):

    if r[0]<7:
        return 'si'
    else:
        return 'air'

    isInsideSol = np.sqrt(r[0]**2+r[1]**2)<=Params.solRad and abs(r[2])<Params.solLength/2
    if isInsideSol:
        mat = 'fe'
    else:
        mat = 'fe'
    return mat
    

def getScatteringParams(x, dt):

    mat = getMaterial(x[:3])

    Z,A,rho,X0 = Params.materials[mat]

    z = Params.Q

    p = x[3:]
    magp = np.linalg.norm(p)
    v = p/np.sqrt(magp**2 + Params.m**2)
    beta = np.linalg.norm(v)

    ds = beta * 2.9979e1 * dt

    Xc = np.sqrt(0.1569 * z**2 * Z*(Z+1) * rho * ds / (magp**2 * beta**2 * A))
    b = np.log(6700*z**2*Z**(1./3)*(Z+1)*rho*ds/A / (beta**2+1.77e-4*z**2*Z**2))

    ## we want to solve the equation B-log(B) = b. Using Newton-Raphson

    B = b
    prevB = 2*B
    
    f = lambda x: x-np.log(x)-b
    fp = lambda x: 1-1./x

    while abs((B-prevB)/prevB)>0.001:
        prevB = B
        B = B - f(B)/fp(B)

        
    # use B+1 for correction at intermediate angles
    return Xc, B+1


def getBField(r):

    if Params.BMag==0:
        return np.zeros(3)
    
    r = np.array(r)
    m = np.array([0.,0.,1500. * Params.BMag/3.856])
    
    isInsideSol = np.sqrt(r[0]**2+r[1]**2)<=Params.solRad and abs(r[2])<Params.solLength/2

    if abs(r[2])<Params.solLength/2:
        if isInsideSol:
            r = np.array([Params.solRad,0,0])
        else:
            r = np.array([np.sqrt(r[0]**2+r[1]**2), 0, 0])
            #r = np.array([Params.solRad,0,0])
    elif r[2]<-Params.solLength/2:
        r[2] = r[2]+Params.solLength/2
    elif r[2]>Params.solLength/2:
        r[2] = r[2]-Params.solLength/2

    mag = np.linalg.norm(r)
    # if mag<Params.solRad:
    #     mag = Params.solRad
        
    B = 3.0*r*np.dot(m,r)/mag**5 - m/mag**3

    if isInsideSol:
        B = -B

    return B

