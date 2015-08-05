#! /usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import Params
import Integrator
import Detector
import Drawing

# VIS or STATS
mode = "VIS"
visWithStats = False

if mode=="VIS":
    ntrajs = 15
    trajs = []
if mode=="STATS":
    ntrajs = 2500
    trajs = []

outname = "data/detectorHits/stats"

Params.BFieldType = 'cms'
Params.MSCtype = 'PDG'
Params.EnergyLossOn = True
Params.Q = 1
Params.m = 105.

Detector.LoadBField("bfield/bfield.pkl")

# make sure numbers are new each run
ROOT.gRandom.SetSeed(0)

rootfile = ROOT.TFile("pdist/test.root")
p_eta_dist = rootfile.Get("samplehisto")

dt = 0.2
nsteps = 400

# set up detector plane

# normal to the plane (x,0,z)
normToDetect = np.array([5,0,15])
# distance from origin to plane
distToDetect = np.linalg.norm(normToDetect)
normToDetect = normToDetect/np.linalg.norm(normToDetect)

#center point
center = normToDetect*distToDetect
# y-axis (in plane)
detVert = np.array([0,1,0])
# another orthogonal vector to norm in plane ((normToDetect, detVert, detOrth) form an ON basis)
detOrth = np.cross(normToDetect,detVert)

detWidth = 4.
detHeight = 4.

detectorDict = {"norm":normToDetect, "dist":distToDetect, "vert":detVert, 
                "orth":detOrth, "width":detWidth, "height":detHeight}

c1 = center + detOrth*detWidth/2 + detVert*detHeight/2
c2 = center + detOrth*detWidth/2 - detVert*detHeight/2
c3 = center - detOrth*detWidth/2 - detVert*detHeight/2
c4 = center - detOrth*detWidth/2 + detVert*detHeight/2

intersects = []

stats = ROOT.TNtuple("stats","stats","p:pT:eta:phi:theta:thOrth:thVert")
ntotaltrajs = 0

txtfile = open(outname+".txt",'w')
txtfile.close()

while len(trajs)<ntrajs:
    magp = ROOT.Double(1e9)
    eta = ROOT.Double(-1)

    etalow = 1.2
    etahigh = 2.5

    while eta<etalow or eta>etahigh:
        p_eta_dist.GetRandom2(magp,eta)

    th = 2*np.arctan(np.exp(-eta))
    phimin, phimax = -1.4,1.4
    phi = np.random.rand() * (phimax-phimin) + phimin
    Params.Q = np.random.randint(2)*2 - 1

    p = 1000*magp * np.array([np.sin(th)*np.cos(phi),np.sin(th)*np.sin(phi),np.cos(th)])
    x0 = np.array([0,0,0,p[0],p[1],p[2]])

    traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
    ntotaltrajs += 1
    if mode=="VIS":
        trajs.append(traj)

    intersection, theta, thOrth, thVert = Detector.FindIntersection(traj, detectorDict)
    if intersection != None:
        intersects.append(intersection)
        print len(trajs), ": p =",magp, ", eta =", eta, ", phi =", phi
        if mode=="VIS":
            pass
        if mode=="STATS":
            if visWithStats:
                trajs.append(traj)
            else:
                trajs.append(0)
            stats.Fill(magp,magp*np.sin(th),eta,phi,theta, thOrth, thVert)
            txtfile = open(outname+".txt",'a')
            txtfile.write("{0:f}\t{1:f}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\n".format(magp,magp*np.sin(th),eta,phi,theta,thOrth,thVert))

fid = ROOT.TFile(outname+".root","RECREATE")
stats.Write()
fid.Close()

print "Efficiency:", float(len(intersects))/ntotaltrajs

if mode=="VIS" or visWithStats:
    plt.figure(num=1, figsize=(15,7))

    Drawing.Draw3Dtrajs(trajs, subplot=121)

    Drawing.DrawLine(c1,c2,is3d=True)
    Drawing.DrawLine(c2,c3,is3d=True)
    Drawing.DrawLine(c3,c4,is3d=True)
    Drawing.DrawLine(c4,c1,is3d=True)

    for i in range(len(intersects)):
        Drawing.DrawLine(intersects[i],intersects[i],is3d=True,linestyle='None',marker='o',color='r')

    Drawing.DrawXYslice(trajs, subplot=122)

    plt.figure(num=2)
    Drawing.DrawXZslice(trajs)

    plt.show()
