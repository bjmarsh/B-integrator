#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#fname = "DY_v3.txt"
#fname = "QCD_v3.txt"
fname = "WJets_v3.txt"

with open(fname) as fid:
    data = np.loadtxt(fid)
data[:,5:8] *= 180/np.pi

pcut = 0.

goodrows = data[:,1]>pcut
data = data[goodrows,:]

print "Number of events:", data.shape[0]

dataP = data[data[:,0]>0, 1:]
dataN = data[data[:,0]<0,1:]

nbins = 50

plt.figure(num=1, figsize=(19,9))

plt.subplot(2,5,1)
plt.hist((dataP[:,0],dataN[:,0]), range=(0,200), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('p')

plt.subplot(2,5,2)
plt.hist((dataP[:,1],dataN[:,1]), range=(0,200), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('pT')

plt.subplot(2,5,3)
plt.hist((dataP[:,2],dataN[:,2]), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('eta')

plt.subplot(2,5,4)
plt.hist((dataP[:,3],dataN[:,3]), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('phi')

plt.subplot(2,5,6)
plt.hist((dataP[:,4],dataN[:,4]), range=(0,15), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('theta')

plt.subplot(2,5,7)
plt.hist((dataP[:,5],dataN[:,5]), range=(-10,10), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('thetaW')

plt.subplot(2,5,8)
plt.hist((dataP[:,6],dataN[:,6]), range=(-10,10), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('thetaV')

plt.subplot(2,5,9)
plt.hist((dataP[:,7],dataN[:,7]), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('w')

plt.subplot(2,5,10)
plt.hist((dataP[:,8],dataN[:,8]), bins=nbins, histtype='step', stacked=False, color=('b','r'))
plt.title('v')


plt.figure(2)

#phi vs. v
maxp = np.amax(data[:,1])
maxp = 100
plt.scatter(data[:,9],data[:,4], s=15, c=data[:,1]*data[:,0], cmap='Spectral', vmin=-maxp, vmax=maxp, alpha=1.0, linewidth=0.3, edgecolor=(0.5,0.5,0.5))
#plt.scatter(dataN[:,8],dataN[:,3], s=15, c=-dataN[:,0], vmax=np.max(dataN[:,0]), cmap='Spectral', alpha=0.5, linewidth=0.3, edgecolor=(0.5,0.5,0.5))
plt.xlabel('v')
plt.ylabel('phi')
# plt.gca().set_xlim(-1,1)
# plt.gca().set_ylim(-0.4,0.4)

# # # thetaV vs pT
# plt.scatter(dataP[:,6],dataP[:,1], s=10, c='b', alpha=0.3, linewidth=0)
# plt.scatter(dataN[:,6],dataN[:,1], s=10, c='r', alpha=0.3, linewidth=0)
# plt.xlabel('thetaV (deg)')
# plt.ylabel('pT (GeV)')
# plt.gca().set_ylim(-10,300)

# # Eta vs pT
# plt.scatter(dataP[:,1],dataP[:,2], c='b', alpha=0.5, linewidth=0)
# plt.scatter(dataN[:,1],dataN[:,2], c='r', alpha=0.5, linewidth=0)

# # Eta vs phi
# plt.scatter(dataP[:,3],dataP[:,2], c='b', alpha=0.5, linewidth=0)
# plt.scatter(dataN[:,3],dataN[:,2], c='r', alpha=0.5, linewidth=0)

# # Eta vs w
# plt.scatter(dataP[:,7],dataP[:,2], c='b', alpha=0.5, linewidth=0)
# plt.scatter(dataN[:,7],dataN[:,2], c='r', alpha=0.5, linewidth=0)
# plt.xlabel('w')
# plt.ylabel('eta')

# # w vs v
# plt.scatter(dataP[:,8],dataP[:,7], s=10, c='b', alpha=0.5, linewidth=0)
# plt.scatter(dataN[:,8],dataN[:,7], s=10, c='r', alpha=0.5, linewidth=0)
# plt.xlabel('v')
# plt.ylabel('w')

# # thetaV vs. v
# plt.scatter(dataP[:,8],dataP[:,6], s=10, c='b', alpha=0.5, linewidth=0)
# plt.scatter(dataN[:,8],dataN[:,6], s=10, c='r', alpha=0.5, linewidth=0)
# plt.xlabel('v')
# plt.ylabel('thetaV')

# # p vs v
# plt.scatter(dataP[:,8],dataP[:,0], s=10, c='b', alpha=0.5, linewidth=0)
# plt.scatter(dataN[:,8],dataN[:,0], s=10, c='r', alpha=0.5, linewidth=0)
# plt.xlabel('v')
# plt.ylabel('p')


thetaW = data[:,6]
mean = np.mean(thetaW)
std = np.std(thetaW)
N = thetaW.size

print "thetaW: mean =", mean, ", std =", std

plt.show()


