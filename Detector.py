## Detector.py
## methods relating to detector and environment properties

import numpy as np
import Params

def LoadBField(fname):
    
    x = np.arange(-Params.RMAX, Params.RMAX+1e-10, Params.DR)
    z = np.arange(Params.ZMIN, Params.ZMAX+1e-10, Params.DZ)
    
    Z,X = np.meshgrid(z,x)

    Params.Bx = np.zeros(X.shape)
    Params.By = np.zeros(X.shape)
    Params.Bz = np.zeros(X.shape)
    Params.Bmag = np.zeros(X.shape)
    
    with open(fname,'r') as fid:
        for line in fid:
            sp = line.strip().split()
            if len(sp)!=4:
                continue

            r = float(sp[0])
            z = float(sp[1])
            phi = float(sp[2])
            B = tuple([float (a) for a in sp[3].strip("()").split(",")])

            if (phi==0 or phi==180):
                x = r
                if phi==180:
                    x*=-1
                iz = int((z+Params.ZMAX)/Params.DZ)
                ix = int((x+Params.RMAX)/Params.DR)

            Params.Bmag[ix,iz] = np.linalg.norm(B)
            Params.Bx[ix,iz] = B[0]
            Params.By[ix,iz] = B[1]
            Params.Bz[ix,iz] = B[2]


def getMaterial(r):

    isInsideSol = np.sqrt(r[0]**2+r[1]**2)<=Params.solRad and abs(r[2])<Params.solLength/2
    if isInsideSol:
        mat = 'fe'
    else:
        mat = 'air'
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


def getBField(x,y,z):

    ## correct for cm usage in bfield file
    x *= 100
    y *= 100
    z *= 100
    r = np.sqrt(x**2+y**2)

    if z>Params.ZMIN and z<Params.ZMAX and r<Params.RMAX:

        r = np.sqrt(x**2+y**2)
        
        nearR = int(Params.DR*round(r/Params.DR))
        nearZ = int(Params.DZ*round(z/Params.DZ))
        
        ir = (nearR+Params.RMAX)/Params.DR
        iz = (nearZ+Params.ZMAX)/Params.DZ
        
        Bx = Params.Bx[ir,iz]
        By = Params.By[ir,iz]
        Bz = Params.Bz[ir,iz]
    
        return np.array([Bx,By,Bz])

    else:
        return np.zeros(3)

    # if Params.BMag==0:
    #     return np.zeros(3)
    
    # r = np.array(r)
    # m = np.array([0.,0.,1500. * Params.BMag/3.856])
    
    # isInsideSol = np.sqrt(r[0]**2+r[1]**2)<=Params.solRad and abs(r[2])<Params.solLength/2

    # if abs(r[2])<Params.solLength/2:
    #     if isInsideSol:
    #         r = np.array([Params.solRad,0,0])
    #     else:
    #         r = np.array([np.sqrt(r[0]**2+r[1]**2), 0, 0])
    #         #r = np.array([Params.solRad,0,0])
    # elif r[2]<-Params.solLength/2:
    #     r[2] = r[2]+Params.solLength/2
    # elif r[2]>Params.solLength/2:
    #     r[2] = r[2]-Params.solLength/2

    # mag = np.linalg.norm(r)
    # # if mag<Params.solRad:
    # #     mag = Params.solRad
        
    # B = 3.0*r*np.dot(m,r)/mag**5 - m/mag**3

    # if isInsideSol:
    #     B = -B

    # return B

