#! /usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import Params
import Integrator
import Detector
import Drawing

Params.BFieldType = 'cms'
Params.MSCtype = 'kuhn'
Params.Q = 1
Params.m = 105.

Detector.LoadBField("bfield/bfield.pkl")

ROOT.gRandom.SetSeed(0)

rootfile = ROOT.TFile("pdist/test.root")
p_eta_dist = rootfile.Get("samplehisto")

dt = 0.2
nsteps = 2000

ntrajs = 5
trajs = [None for i in range(ntrajs)]

for i in range(ntrajs):
    magp = ROOT.Double()
    eta = ROOT.Double()

    p_eta_dist.GetRandom2(magp,eta)

    print magp, eta

    th = 2*np.arctan(np.exp(-eta))
    phi = np.random.rand() * 2*np.pi

    p = 1000*magp * np.array([np.sin(th)*np.cos(phi),np.sin(th)*np.sin(phi),np.cos(th)])
    x0 = np.array([0,0,0,p[0],p[1],p[2]])

    trajs[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)


plt.figure(num=1, figsize=(15,7))

Drawing.Draw3Dtrajs(trajs, subplot=121)
Drawing.DrawXYslice(trajs, subplot=122)

plt.show()
