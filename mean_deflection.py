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

# modify globals from their default values
Params.BMag = 0.0

dt = 0.2
nsteps = 1000

Params.MSCtype = 'Kuhn'

pvals = np.arange(1000,20001,1000)

meanrvals = []
sdevrvals = []
meanyvals = []
sdevyvals = []

Nsamp = 100

for p in pvals:

    p0 = [p,0,0]
    x0 = np.array([0,0,0]+p0)
    
    yvals = []
    rvals = []
    thvals = []
    for i in range(Nsamp):
        traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, cutoff=5, cutoffaxis=0)
        
        # get the coordinates as it leaves material
        k = (5-traj[0,-2])/(traj[0,-1]-traj[0,-2])
        cross = traj[:3,-2]+k*(traj[:3,-1]-traj[:3,-2])

        rvals.append(np.sqrt(cross[1]**2 + cross[2]**2))
        yvals.append(cross[1])
        
        # get the projected angle as it exits material
        projth = np.arctan2(traj[1,-1],traj[0,-1])
        thvals.append(projth)

    fid = open("data/pt"+str(p),"w")
    for i in range(Nsamp):
        fid.write("{0}\t{1}\t{2}\n".format(rvals[i],yvals[i],thvals[i]))
    fid.close()

    meanrvals.append(np.mean(rvals))
    sdevrvals.append(np.std(rvals))
    meanyvals.append(np.mean(np.absolute(yvals)))
    sdevyvals.append(np.mean(np.std(yvals)))

    if p==5000:
        testth = thvals

meanrvals = np.array(meanrvals)
sdevrvals = np.array(sdevrvals)
meanyvals = np.array(meanyvals)
sdevyvals = np.array(sdevyvals)

fid = open("data/mean_defl","w")
for i in range(len(pvals)):
    fid.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(pvals[i],meanrvals[i],sdevrvals[i],meanyvals[i],sdevyvals[i]))
fid.close()

#plt.errorbar(pvals/1000, meanrvals, yerr=sdevrvals/np.sqrt(Nsamp), fmt='ok')
#plt.errorbar(pvals/1000, meanyvals, yerr=sdevyvals/np.sqrt(Nsamp), fmt='or')

plt.hist(testth, bins=20)        

plt.show()
