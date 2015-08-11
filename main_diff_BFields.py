#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import Integrator
import Detector
import Params
import Drawing

Detector.LoadCoarseBField("bfield/bfield_coarse.pkl")
print "Loaded coarse field"
Detector.LoadFineBField("bfield/bfield_x.pkl","bfield/bfield_y.pkl","bfield/bfield_z.pkl")
print "Loaded fine field"

Params.BFieldType = 'cms'
Params.Q = 1.0
Params.MSCtype = 'none'
Params.EnergyLossOn=False

dt = 0.2
nsteps = 200

p0 = [10000, 0, 0]
x0 = np.array([0,0,0]+p0)

UseFineBField = False
