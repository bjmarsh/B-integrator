#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import Integrator
import Detector
import Params
import Drawing

#Detector.LoadBField("bfield/bfield.pkl")

Params.BFieldType = 'none'
Params.Q = 1.0
scatterType = 'PDG'
Params.MSCtype = scatterType
Params.MatSetup = 'iron'
Params.EnergyLossOn=False

#print "\nDone loading B-field...\nB-field at center:", Detector.getBField(0,0,0), '\n'


# define the initial momenta (in MeV)
init_p = []
#init_p.append([1000, 0, 0])
# init_p.append([4800, 0, 2000])
# init_p.append([10000, 0, 0])
# init_p.append([10000, 5000, 17000])
# init_p.append([30000, 40000, 5000])

init_p.append([2000,0,0])
init_p.append([3000,0,0])
init_p.append([5000,0,0])
init_p.append([10000,0,0])
init_p.append([20000,0,0])

colors = ['r', 'g', 'b', 'c', 'm']

print 'Initial Momenta (colors r,g,b,c):'
for i in range(len(init_p)):
    print ' -',round(np.linalg.norm(init_p[i])/1000, 1), "GeV"

dt = 0.1
nsteps = 700

trajs = list(range(len(init_p)))
trajs_noMSC = list(range(len(init_p)))

for i in range(len(init_p)):
    x0 = np.array([0,0,0]+init_p[i])
    Params.MSCtype = 'none'
    trajs_noMSC[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)
    Params.MSCtype = scatterType
    trajs[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)

plt.figure(1)

time = np.arange(0,dt*nsteps+1e-10, dt)
for i in range(len(init_p)):
    maxI = min(trajs[i].shape[1], trajs_noMSC[i].shape[1])
    ydiff = trajs[i][1,:maxI]-trajs_noMSC[i][1,:maxI]
    zdiff = trajs[i][2,:maxI]-trajs_noMSC[i][2,:maxI]
    plt.plot(trajs[i][0,:maxI], zdiff, color=colors[i], 
             label=str(int(init_p[i][0]/1000))+" GeV")

plt.xlabel('Horizontal distance (m)')
plt.ylabel('Displacement (m)')
plt.title('Paths through block of iron')
plt.legend(loc='upper right',fontsize='small')

Params.MatSetup = 'sife'
Params.EnergyLossOn = True

for i in range(len(init_p)):
    x0 = np.array([0,0,0]+init_p[i])
    trajs[i] = Integrator.rk4(x0, Integrator.traverseBField, dt, nsteps)

plt.figure(2)

for i in range(len(init_p)):
    trajE = np.linalg.norm(trajs[i][3:,:],axis=0)/1000.
    plt.plot(trajs[i][0,:],trajE, color=colors[i],
             label=str(int(init_p[i][0]/1000))+" GeV")

plt.plot([4,4],[0,20],'--k')
plt.text(3,13,'Si',size='large')
plt.text(4.6,13,'Fe',size='large')

plt.legend()
plt.xlabel('Horizontal Distance (m)')
plt.ylabel('Energy (GeV)')
plt.title('Energy vs. distance')
plt.show()
exit(0)

fig = plt.figure(2, figsize=(14.5,5.5))
## xz slice
x = np.arange(-Params.RMAX, Params.RMAX+1e-10, Params.DR)/100
z = np.arange(Params.ZMIN, Params.ZMAX+1e-10, Params.DZ)/100
    
Z,X = np.meshgrid(z,x)

plt.subplot2grid((1,5),(0,0),colspan=3)

# draw mag field

if Params.BFieldLoaded:
    mag = np.append(Params.Bmag[::-1,:,0],Params.Bmag[1:,:,180/Params.DPHI],0)
    bmplot = plt.pcolor(Z,X,mag,cmap='afmhot',vmax=5.0)
    #bmcb = plt.colorbar(bmplot, orientation='horizontal')
    
    #if np.linalg.norm(Detector.getBField(0,0,0)) != 0:
    #    plt.quiver(Z,X,Bzy,Bxy)

sl = Params.solLength
sr = Params.solRad

# draw trajectory
for i in range(len(init_p)):
    plt.plot(trajs[i][2,:],trajs[i][0,:],'-', linewidth=2, color=colors[i])

plt.axis([-15,15,-9,9])
plt.xlabel("z (m)")
plt.ylabel("x (m)")

## xy slice

plt.subplot2grid((1,5),(0,3), colspan=2)

Drawing.DrawXYslice(trajs,ax=plt.gca())

# # draw dotted lines to make figure for slides
# pf = trajs[0][3:,-1]
# xf = trajs[0][:3,-1]
# t = -np.dot(pf,xf)/np.dot(pf,pf)
# xi = xf+1.4*t*pf
# plt.plot([xi[0],xf[0]],[xi[1],xf[1]],'--g')
# plt.plot([-4,9],[0,0],'--g')
# xi = xf+t*pf
# plt.plot([0,xi[0]],[0,xi[1]],'-m')
# plt.text(-0.08,-0.6,r'b', size='large')

plt.axis([-9,9,-9,9])
plt.xlabel("x (m)")
plt.ylabel("y (m)")

## 3d view

fig = plt.figure(num=3, figsize=(8, 8))

Drawing.Draw3Dtrajs(trajs)


plt.show()



