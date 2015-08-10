## Detector.py
## methods relating to detector and environment properties

import cPickle as pickle
import numpy as np
import Params

def LoadBField(fname):
    
    if Params.BFieldUsePickle:
        Params.Bx,Params.By,Params.Bz,Params.Bmag = pickle.load(open(fname,"rb"))
        Params.BFieldLoaded = True
        return

    NR = Params.RMAX/Params.DR + 1
    NZ = (Params.ZMAX-Params.ZMIN)/Params.DZ + 1
    NPHI = (Params.PHIMAX-Params.PHIMIN)/Params.DPHI + 1


    Params.Bx   = np.zeros((NR,NZ,NPHI))
    Params.By   = np.zeros((NR,NZ,NPHI))
    Params.Bz   = np.zeros((NR,NZ,NPHI))
    Params.Bmag = np.zeros((NR,NZ,NPHI))
    
    with open(fname,'r') as fid:
        for line in fid:
            sp = line.strip().split()
            if len(sp)!=4:
                continue

            r = float(sp[0])
            z = float(sp[1])
            phi = float(sp[2])
            B = tuple([float (a) for a in sp[3].strip("()").split(",")])

            iz = int((z-Params.ZMIN)/Params.DZ)
            ir = int(r/Params.DR)
            iphi = int((phi-Params.PHIMIN)/Params.DPHI)

            Params.Bmag[ir,iz,iphi] = np.linalg.norm(B)
            Params.Bx[ir,iz,iphi] = B[0]
            Params.By[ir,iz,iphi] = B[1]
            Params.Bz[ir,iz,iphi] = B[2]

    Params.BFieldLoaded = True



def getMaterial(x,y,z):

    if Params.MatSetup == 'iron':
        return 'fe'

    if Params.MatSetup == 'sife':
        if x<4:
            return 'si'
        else:
            return 'fe'

    withinLength = -Params.solLength/2 < z < Params.solLength/2
    r = np.sqrt(x**2+y**2)

    if not withinLength:
        return 'air'
    
    if r < 1.3:
        mat = 'si'
    elif r < 1.8:
        mat = 'pbwo4'
    elif r < 2.95:
        mat = 'fe'
    elif r < 4.0:
        mat = 'fe'
    elif r < 7.0:
        mat = 'fe'
    else:
        mat = 'air'

    return mat
    

def getBField(x,y,z):

    if Params.BFieldType.lower() == 'none':
        return np.zeros(3)

    if Params.BFieldType.lower() == 'uniform':
        return np.array([0.,0.,1.])

    if Params.BFieldType.lower() == 'updown':
        r = np.sqrt(x**2+y**2)
        if r<4:
            return np.array([0,0,3])
        else:
            return np.array([0,0,-3])/r

    ## correct for cm usage in bfield file
    x *= 100
    y *= 100
    z *= 100
    r = np.sqrt(x**2+y**2)

    if z>Params.ZMIN and z<Params.ZMAX and r<Params.RMAX:

        r = np.sqrt(x**2+y**2)
        phi = np.arctan2(y,x) * 180/np.pi

        if phi<0:
            phi += 360
        
        nearR = int(Params.DR*round(r/Params.DR))
        nearZ = int(Params.DZ*round(z/Params.DZ))
        nearPHI = int(Params.DPHI*round(phi/Params.DPHI))
        
        if nearPHI==360:
            nearPHI = 0

        ir = nearR/Params.DR
        iz = (nearZ-Params.ZMIN)/Params.DZ
        iphi = (nearPHI-Params.PHIMIN)/Params.DPHI
        
        Bx = Params.Bx[ir,iz,iphi]
        By = Params.By[ir,iz,iphi]
        Bz = Params.Bz[ir,iz,iphi]
    
        return np.array([Bx,By,Bz])

    else:
        return np.zeros(3)

def FindIntersection(traj, detectorDict):
    # find the intersection with a plane with normal norm
    # and distance to origin dist. returns None if no intersection

    norm = detectorDict["norm"]
    dist = detectorDict["dist"]

    for i in range(traj.shape[1]-1):
        p1 = traj[:3,i]
        p2 = traj[:3,i+1]
        
        proj2 = np.dot(p2,norm)

        if proj2>=dist:
            proj1 = np.dot(p1,norm)
            intersect = p1+(dist-proj1)/(proj2-proj1)*(p2-p1)

            vHat = detectorDict["v"]
            wHat = detectorDict["w"]
            center = norm*dist

            w = np.dot(intersect,wHat)
            v = np.dot(intersect,vHat)
            
            if abs(w) < detectorDict["width"]/2 and \
               abs(v) < detectorDict["height"]/2:
                unit = (p2-p1)/np.linalg.norm(p2-p1)
                theta = np.arccos(np.dot(unit,norm))
                
                projW = np.dot(unit,wHat)
                projV = np.dot(unit,vHat)

                thW = np.arcsin(projW/np.linalg.norm(unit-projV*vHat))
                thV = np.arcsin(projV/np.linalg.norm(unit-projW*wHat))

                return intersect,theta,thW,thV

            break

    return None, None, None, None












