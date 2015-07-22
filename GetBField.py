import numpy as np
from Params import *

def getBField(r):

    global solRad, solLength

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
