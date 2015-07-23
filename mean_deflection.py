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
    for i in range(Nsamp):
        traj = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps, cutoff=7, cutoffaxis=0)
        rvals.append(np.sqrt(traj[1,-1]**2 + traj[2,-1]**2))
        yvals.append(traj[1,-1])
        
    yfid = open("data/y_pt"+str(p),"w")
    yfid.write('\n'.join([str(a) for a in yvals]))
    yfid.close()
    rfid = open("data/r_pt"+str(p),"w")
    rfid.write('\n'.join([str(a) for a in rvals]))
    rfid.close()

    meanrvals.append(np.mean(rvals))
    sdevrvals.append(np.std(rvals))
    meanyvals.append(np.mean(np.absolute(yvals)))
    sdevyvals.append(np.mean(np.std(yvals)))

meanrvals = np.array(meanrvals)
sdevrvals = np.array(sdevrvals)
meanyvals = np.array(meanyvals)
sdevyvals = np.array(sdevyvals)

fid = open("data/mean_defl","w")
for i in range(len(pvals)):
    fid.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(pvals[i],meanrvals[i],sdevrvals[i],meanyvals[i],sdevyvals[i]))
fid.close()

plt.errorbar(pvals/1000, meanrvals, yerr=sdevrvals/np.sqrt(Nsamp), fmt='ok')
plt.errorbar(pvals/1000, meanyvals, yerr=sdevyvals/np.sqrt(Nsamp), fmt='or')
        
# plt.hist(yvals, bins=20)

plt.show()
