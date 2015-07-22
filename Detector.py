import numpy as np
from Params import *

def getScatteringParams(x, dt):
    return


def getRadiationLength(r):
    # return the radiation length of the material at position r
    # for now, assumes constant material within solenoid, vacuum outside

    global solRad, solLength

    isInsideSol = np.sqrt(r[0]**2+r[1]**2)<=solRad and abs(r[2])<solLength/2
    
    if isInsideSol:
        return 0.05
    else:
        return 0

def getBField(r):

    global solRad, solLength, BMag

    r = np.array(r)
    m = np.array([0.,0.,1500. * BMag/3.856])
    
    isInsideSol = np.sqrt(r[0]**2+r[1]**2)<=solRad and abs(r[2])<solLength/2

    if abs(r[2])<solLength/2:
        if isInsideSol:
            r = np.array([solRad,0,0])
        else:
            r = np.array([np.sqrt(r[0]**2+r[1]**2), 0, 0])
            #r = np.array([solRad,0,0])
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

