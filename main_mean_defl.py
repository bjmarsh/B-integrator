#! /usr/bin/python

import timeit
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import Integrator
import Detector
import MultipleScatter
import Params

Detector.LoadBField("bfield/bfield.txt")

dt = 0.2
nsteps = 1000

pvals = np.arange(2000,20001,1000)

meanrvals = []
sdevrvals = []
meanyvals = []
sdevyvals = []

Nsamp = 10

#function to get r value of a point
getr = lambda x: np.sqrt(x[0]**2+x[1]**2)

for p in pvals:

    print p

    p0 = [p,0,0]
    x0 = np.array([0,0,0]+p0)
    
    rvals = []
    thvals = []

    ## get the trajectory & end coords with no multiple scattering
    Params.MSCtype = 'none'
    traj_noMSC = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, cutoff=8, cutoffaxis=3)
    k = (8-getr(traj_noMSC[:,-2]))/(getr(traj_noMSC[:,-1])-getr(traj_noMSC[:,-2]))
    cross_noMSC = traj_noMSC[:3,-2]+k*(traj_noMSC[:3,-1]-traj_noMSC[:3,-2])
    projth_noMSC = np.arctan2(traj_noMSC[1,-1],traj_noMSC[0,-1])
    
    Params.MSCtype = 'Kuhn'

    for i in range(Nsamp):

        traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, cutoff=8, cutoffaxis=3)
        
        # get the coordinates as it leaves material
        k = (8-getr(traj[:,-2]))/(getr(traj[:,-1])-getr(traj[:,-2]))
        cross = traj[:3,-2]+k*(traj[:3,-1]-traj[:3,-2])

        rvals.append(np.sqrt((cross[1]-cross_noMSC[1])**2 + (cross[2]-cross_noMSC[2])**2))
        
        # get the projected angle as it exits material
        projth = np.arctan2(traj[1,-1],traj[0,-1])
        thvals.append(abs(projth-projth_noMSC))

    fid = open("data/p"+str(p),"w")
    for i in range(Nsamp):
        fid.write("{0}\t{1}\n".format(rvals[i],thvals[i]))
    fid.close()

    meanrvals.append(np.mean(rvals))
    sdevrvals.append(np.std(rvals))


meanrvals = np.array(meanrvals)
sdevrvals = np.array(sdevrvals)

fid = open("data/mean_defl","w")
for i in range(len(pvals)):
    fid.write("{0}\t{1}\t{2}\n".format(pvals[i],meanrvals[i],sdevrvals[i]))
fid.close()

#plt.errorbar(pvals/1000, meanrvals, yerr=sdevrvals/np.sqrt(Nsamp), fmt='ok')
#plt.errorbar(pvals/1000, meanyvals, yerr=sdevyvals/np.sqrt(Nsamp), fmt='or')

#plt.show()
