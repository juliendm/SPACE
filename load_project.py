#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
from optparse import OptionParser
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors


# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    ndim = 11
    lift_model = Surfpack('LIFT',ndim)

    project = SPACE.io.load_data('RESPONSE_SURFACE_DV_SUP/project.pkl')
    project.compile_designs()

    for index in range(len(project.designs)):
 
        design = project.designs[index].design
    
        funcs = design.funcs
        config = design.config

        dv_mach = float(config.MACH_NUMBER)
        dv_re = np.log10(float(config.REYNOLDS_NUMBER)) + 3.0/8.0*dv_mach - 7.0
        if dv_mach >= 1.0:
            dv_aoa = (float(config.AoA) - (3.92857*dv_mach+3.57143)) / (1.78571*dv_mach+5.71429) # Supersonic
        else:
            dv_aoa = float(config.AoA)/7.5-1.0 # Subsobic
        dv_bf = float(config.BODY_FLAP_DEF)*np.pi/180.0
        dv_el = float(config.ELEVON_DEF)*np.pi/180.0
        dv_geo1 = float(config.DV1)
        # # dv2 = dv2_real + dv4
        dv_geo2 = float(config.DV2) - float(config.DV4) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dv_geo3 = float(config.DV3)
        dv_geo4 = float(config.DV4)
        dv_geo5 = float(config.DV5)
        dv_geo6 = float(config.DV6)

        dvs = [dv_mach,dv_re,dv_aoa,dv_bf,dv_el,dv_geo1,dv_geo2,dv_geo3,dv_geo4,dv_geo5,dv_geo6]

        if hasattr(funcs,'LIFT'):
            lift_model.add(dvs,funcs.LIFT)


    lift_model.save_data('build_points_lift.dat')

    lift_model.build('kriging')
    lift_model.save_model('lift.sps')

    # lift_model.load_model('lift.sps')

    nx = 50

    fig = plt.figure()
    ax = Axes3D(fig)


    model = lift_model
    model_range = [-0.1,0.6]

    mach_range = [[ 1.1, 1.5],[ 1.5, 2.0],[ 2.0, 2.5],[ 2.5, 3.0],[ 3.0, 4.0],[ 4.0, 5.0],[ 5.0, 8.0]]
    aoa_range  = [[ 0.0,15.0],[ 2.5,20.0],[ 5.0,25.0],[ 7.5,27.5],[10.0,30.0],[12.5,40.0],[15.0,50.0]]

    for index in range(len(mach_range)):

        MACH = np.linspace(mach_range[index][0], mach_range[index][1], nx)
        AOA = np.linspace(aoa_range[index][0], aoa_range[index][1], nx)
        MACH,AOA = np.meshgrid(MACH,AOA)
        COEFF = MACH*0.

        for ix in range(nx):
            for iy in range(nx):

                dv_mach = MACH[ix,iy]

                if dv_mach >= 1.0:
                    dv_aoa = (AOA[ix,iy] - (3.92857*dv_mach+3.57143)) / (1.78571*dv_mach+5.71429) # Supersonic
                else:
                    dv_aoa = AOA[ix,iy]/7.5-1.0 # Subsobic

                COEFF[ix,iy] = model.eval([dv_mach,0.0,dv_aoa,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

        ax.plot_surface(MACH,AOA,COEFF, norm=colors.Normalize(vmin=model_range[0],vmax=model_range[1]), cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)


    # for ix in range(nx):
    #     for iy in range(nx):
    #         ZP[ix,iy] = lift_model.eval([XP[ix,iy],1E6,YP[ix,iy],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    # ax.plot_surface(XP,YP,ZP, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)

    # for ix in range(nx):
    #     for iy in range(nx):
    #         ZP[ix,iy] = lift_model.eval([XP[ix,iy],1E5,YP[ix,iy],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    # ax.plot_surface(XP,YP,ZP, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)

    # for ix in range(nx):
    #     for iy in range(nx):
    #         ZP[ix,iy] = lift_model.eval([XP[ix,iy],4E6,YP[ix,iy],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    # ax.plot_surface(XP,YP,ZP, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)

    # for ix in range(nx):
    #     for iy in range(nx):
    #         ZP[ix,iy] = lift_model.eval([XP[ix,iy],5E6,YP[ix,iy],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    # ax.plot_surface(XP,YP,ZP, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)

    plt.xlabel('Mach', fontsize=18)
    plt.ylabel('AoA', fontsize=18)

    fig.savefig('fig.png')   # save the figure to file
    plt.close(fig)    # close the figure

    # plt.draw()
    # plt.show()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
