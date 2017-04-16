#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np

from Surfpack.surfpack import Surfpack

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

def main():

    objective = 'lift'
#    objective = 'drag'

    nx = 80

    XB_sub = np.array([[0.5, 0.9],
       [-5.0, 15.0],
       [-0.5, 0.5],
       [-0.5, 0.5],
       [-0.5, 0.5]])

    XB_sup = np.array([[1.2, 9.0],
       [0.0, 20.0],
       [-0.5, 0.5],
       [-0.5, 0.5],
       [-0.5, 0.5]])

    ndim = len(XB_sup)

#    surfpack = Surfpack(ndim)
#    surfpack.load(objective + '_sub','sps/sp_gp_model.' + objective + '_sub.sps')
#    surfpack.load(objective + '_sup','sps/sp_gp_model.' + objective + '.sps')

    fig = plt.figure()
    ax = Axes3D(fig)

#    XP = np.linspace(XB_sub[0][0], XB_sub[0][1], nx)
#    YP = np.linspace(XB_sub[1][0], XB_sub[1][1], nx)
#    XP,YP = np.meshgrid(XP,YP)
#    ZP = YP*0.
#
#    for ix in range(nx):
#        for iy in range(nx):
#            ZP[ix,iy] = surfpack.eval(objective + '_sub',[XP[ix,iy],YP[ix,iy],-0.5,0.5,-0.5])
#    surf = ax.plot_surface(XP,YP,ZP, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)


#    XP = np.linspace(XB_sup[0][0], XB_sup[0][1], nx)
#    YP = np.linspace(XB_sup[1][0], XB_sup[1][1], nx)
#    XP,YP = np.meshgrid(XP,YP)
#    ZP = YP*0.
#
#    for ix in range(nx):
#        for iy in range(nx):
#            ZP[ix,iy] = surfpack.eval(objective + '_sup',[XP[ix,iy],YP[ix,iy],0.0,0.0,0.0])
#    surf = ax.plot_surface(XP,YP,ZP, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)





    surfpack = Surfpack(2)

    XP = np.linspace(-0.5, 0.5, nx)
    YP = np.linspace(-0.5, 0.5, nx)
    XP,YP = np.meshgrid(XP,YP)
    ZP = YP*0.

    data = 'perfo'
    surfpack.load(data,'sps/sp_gp_model.' + data + '.sps')
    for ix in range(nx):
        for iy in range(nx):
            ZP[ix,iy] = surfpack.eval(data,[XP[ix,iy],YP[ix,iy]])
    ax.plot_surface(XP,YP,ZP, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)
    X = np.loadtxt('build_points/enrich_points_' + data + '.dat')
    ax.plot(X[:,0],X[:,1],X[:,2],'k+',mew=3,ms=15)

#    data = 'perfo_one_flylaw'
#    surfpack.load(data,'sps/sp_gp_model.' + data + '.sps')
#    for ix in range(nx):
#        for iy in range(nx):
#            ZP[ix,iy] = surfpack.eval(data,[XP[ix,iy],YP[ix,iy]])
#    ax.plot_surface(XP,YP,ZP, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)
#    X = np.loadtxt('build_points/build_points_' + data + '.dat')
#    ax.plot(X[:,0],X[:,1],X[:,2],'k+',mew=3,ms=15)

    plt.xlabel('dv1', fontsize=18)
    plt.ylabel('dv4', fontsize=18)

    plt.draw()
    plt.show()

#: main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
