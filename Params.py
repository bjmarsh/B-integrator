# Params.py
# control all global variables

import numpy as np

## materials definition. (atomic num, atomic weight, density (g/cm3), radiation length (m))
materials = { "fe"   : (26, 55.845, 7.874, .01757), 
              "si"   : (14, 28.0855, 2.329, .0937),
              "air"  : (7.34, 14.719, 1.205e-3, 3.04e2),  
              "pbwo4": (31.3, 75.8, 8.3, 0.008903)  }

## parameters for dEdx (I, a, k, x0, x1, Cbar, delta0)
dEdx_params = { "fe"   : (286.0, 0.14680, 2.9632, -0.0012, 3.1531, 4.2911, 0.12),
                "si"   : (173.0, 0.14921, 3.2546, 0.2015, 2.8716, 4.4355, 0.14),
                "air"  : (85.7, 0.10914, 3.3994, 1.7418, 4.2759, 10.5961, 0.0),
                "pbwo4": (600.7, 0.22758, 3.0, 0.4068, 3.0023, 5.8528, 0.0) }

## these may be updated in main program
Q = 1  ## in units of e
m = 105.658  ## in MeV 
solRad = 3.6  ## in m
solLength = 15.0   ## in m
MSCtype = 'kuhn'
EnergyLossOn = False
SuppressStoppedWarning = False
BFieldType = 'CMS'
BFieldUsePickle = True
MatSetup = 'cms'

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
