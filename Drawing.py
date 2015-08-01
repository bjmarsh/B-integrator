## Drawing.py
## contains various routines to visualize trajectories

import numpy as np
import matplotlib.pyplot as plt
import Params

## 3d view

def Draw3Dtrajs(trajs, colors=None, ax = None, fig=None, subplot=111):
    # trajs is a list of trajectory arrays as returned by the Integrator.rk4 routine
    # colors is an optional array of colors to use
    # ax is the mpl axes to plot to. If None, createsnew axes
    # fig is the figure to use (if ax not specified). Defaults to plt.gcf()
    # subplot is the subplot to plot to

    from mpl_toolkits.mplot3d import Axes3D

    if ax==None:
        if fig==None:
            fig = plt.gcf()
        ax = fig.add_subplot(subplot, projection='3d')

    if colors==None:
        colors = ['r','g','b','c','m','y']

    nc = len(colors)
        
    for i in range(len(trajs)):
        ax.plot3D(xs=trajs[i][0,:], ys=trajs[i][1,:], zs=trajs[i][2,:], color=colors[i%nc])

    sr, sl = 4., 26.
    t = np.linspace(0, 2*np.pi, 100)
    ax.plot(xs=sr*np.cos(t), ys=sr*np.sin(t), zs=sl/2, color='k')
    ax.plot(xs=sr*np.cos(t), ys=sr*np.sin(t), zs=-sl/2, color='k')
    for i in range(8):
        th = i * 2*np.pi/8
        x = sr*np.cos(th)
        y = sr*np.sin(th)
        ax.plot(xs=[x,x], ys=[y,y], zs=[-sl/2, sl/2], color='k')

    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')

    ax.set_xlim((-9,9))
    ax.set_ylim((-9,9))
    ax.set_zlim((-15,15))



def DrawXYslice(trajs, colors=None, ax=None, fig=None, subplot=111):
    # draws projected trajectories in z=0 plane
    # see above for argument descriptions

    if ax == None:
        if fig==None:
            fig = plt.gcf()
        ax = fig.add_subplot(subplot)

    if colors==None:
        colors = ['r','g','b','c','m','y']

    nc = len(colors)

    # draw solenoid outline
    sr = Params.solRad
    sl = Params.solLength
    t = np.linspace(0, 2*np.pi, 100)
    ax.plot(sr*np.cos(t),sr*np.sin(t), '-')

    # draw trajectory
    for i in range(len(trajs)):
        ax.plot(trajs[i][0,:],trajs[i][1,:],'-', linewidth=2, color=colors[i%nc])

    ax.set_xlim((-9,9))
    ax.set_ylim((-9,9))



