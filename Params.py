# Params.py
# control all global variables

import numpy as np

## materials definition. (atomic num, atomic weight, density (g/cm3), radiation length (m), mean excitation energy (eV))
materials = { "fe"   : (26, 55.845, 7.874, .01757, 173.0), 
              "si"   : (14, 28.0855, 2.329, .0937, 286.0),
              "air"  : (7.34, 14.719, 1.205e-3, 3.04e2, 85.7),  
              "pbwo4": (31.3, 75.8, 8.3, 0.008903, 600.7)  }

## these may be updated in main program
Q = 1  ## in units of e
m = 105.658  ## in MeV 
solRad = 3.5  ## in m
solLength = 21.6   ## in m
MSCtype = 'kuhn'
EnergyLossOn = False
BFieldType = 'CMS'
BFieldUsePickle = True


## internal parameters. don't touch
MSCWarning = False
BFieldLoaded = False

## parameters to load bfield
ZMIN = -1500
ZMAX = 1500
DZ = 10
RMIN = 0
RMAX = 900
DR = 10
PHIMIN = 0
PHIMAX = 355
DPHI = 5
Bx = np.array([])
By = np.array([])
Bz = np.array([])
Bmag = np.array([])
