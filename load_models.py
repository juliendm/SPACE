#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
import scipy as sp
from optparse import OptionParser
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from SPACE.util import DesignVariables
from SPACE.eval import model as spacemodel

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--project", dest="project_folder",
                      help="project folder", metavar="PROJECT_FOLDER")
    parser.add_option("-r", "--regime", dest="regime", default="SUP",
                      help="regime", metavar="REGIME")

    (options, args)=parser.parse_args()

    load_models( options.filename       ,
                 options.project_folder ,
                 options.regime )

#: main()

def load_models( filename         ,
                 project_folder   ,
                 regime = 'SUP' ):

    # Project

    if os.path.exists(project_folder):
        project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
        project.compile_designs()
        config = project.config
    else:
        config = SPACE.io.Config(filename)
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    
    # Start new model

    konfig = copy.deepcopy(config)
    aero = spacemodel.AeroModel(konfig)

    # print aero.static_margin([8.000000e+00,   1.000000e+00,  -1.000000e+00,   2.059314e-01,   5.000000e-01,   5.000000e-01,  -5.000000e-01,   5.000000e-01,  -2.000000e-01,  -5.000000e-01,  -2.000000e-01],11.0)
    # print aero.moment_y([5.710816e+00,   1.000000e+00,  -1.000000e+00,  -2.500000e-01,  -3.500000e-01,  -5.000000e-01,  -5.000000e-01,  -5.000000e-01,  -2.000000e-01,  -5.000000e-01,  -2.000000e-01],11.0)
    # print aero.moment_y([3.665795e+00,   1.000000e+00,  -1.000000e+00,   5.000000e-01,   5.000000e-01,  -5.000000e-01,  -5.000000e-01,   5.000000e-01,  -2.000000e-01,  -5.000000e-01,  -2.000000e-01],11.0)

    # Plot

    nx = 50

    model_range = [-10.0,60.0]

    mach_range = [[ 0.0, 8.0]]
    aoa_range  = [[ 0.0,55.0]]

    cogs = [8.5,9.0,9.5,10.0,10.5,11.0,11.5] # 9.6,9.7,9.8,9.9,

    for cog in cogs:

        print cog

        fig = plt.figure()

        for index in range(len(mach_range)):

            MACH = np.linspace(mach_range[index][0], mach_range[index][1], nx)
            AOA = np.linspace(aoa_range[index][0], aoa_range[index][1], nx)
            MACH, AOA = np.meshgrid(MACH, AOA)
            COEFF = MACH*0.0

            for ix in range(nx):
                for iy in range(nx):
                    dvs = aero.desvar.reverse_variable_change([MACH[ix,iy],1E7,AOA[ix,iy], 0.0,0.0, -0.5,0.5,0.5,0.5, 0.5,0.5])  ##### MUST DO BETTER for Reynolds... Must be the one from trajectory!!! FOR REYNOLDS

                    #COEFF[ix,iy] = aero.trim(dvs,cog)[aero.desvar.trim_index]*180.0/np.pi
                    #COEFF[ix,iy] = aero.static_margin(aero.trim(dvs,cog),cog)
                    COEFF[ix,iy] = aero.max_trimmed_efficiency(dvs,cog)

            #levels = [-float("inf"), -10.0, 30.0, float("inf")]; colors = ('r', 'g', 'r')
            #levels = [-float("inf"), -0.04, float("inf")]; colors = ('r', 'g')
            levels = [-float("inf"), 5.0, float("inf")]; colors = ('r', 'g')

            CS = plt.contourf(MACH, AOA, COEFF, levels, colors=colors)
            plt.clabel(CS, inline=1, fontsize=10)

        plt.xlabel('Mach', fontsize=18)
        plt.ylabel('AoA', fontsize=18)

        plt.ylim([0,55])

#        plt.plot([0,1], [15,15], color='k')
        plt.fill_between([0,1], [15,15], [55,55], color='grey', alpha='0.5', zorder=3)
#        plt.plot([1,8], [15,55], color='k')
        plt.fill_between([1,8], [15,55], [55,55], color='grey', alpha='0.5', zorder=3)

#        plt.plot([1,8], [0 ,15], color='k')
        plt.fill_between([1,8], [0,0], [0,15], color='grey', alpha='0.5', zorder=3)

        fig.savefig(os.path.join(project_folder,'fig_' + str(cog) + '.png'))
        plt.close(fig)


    # # Plot

    # nx = 50

    # fig = plt.figure()
    # ax = Axes3D(fig)

    # model_range = [0.0,0.6]

    # mach_range = [[ 1.1, 1.5],[ 1.5, 2.0],[ 2.0, 2.5],[ 2.5, 3.0],[ 3.0, 4.0],[ 4.0, 5.0],[ 5.0, 8.0]]
    # aoa_range  = [[ 0.0,15.0],[ 2.5,20.0],[ 5.0,25.0],[ 7.5,27.5],[10.0,30.0],[12.5,40.0],[15.0,50.0]]


    # # mach_range = [[ 1.1, 8.0]]
    # # aoa_range  = [[ -5.0,30.0]]

    # for index in range(len(mach_range)):

    #     MACH = np.linspace(mach_range[index][0], mach_range[index][1], nx)
    #     AOA = np.linspace(aoa_range[index][0], aoa_range[index][1], nx)
    #     MACH, AOA = np.meshgrid(MACH, AOA)
    #     COEFF = MACH*0.0

    #     for ix in range(nx):
    #         for iy in range(nx):

    #             dv_mach = MACH[ix,iy]
    #             dv_re, dv_aoa = aero.desvar_sup.reverse_variable_change(dv_mach, 1E7, AOA[ix,iy])
    #             COEFF[ix,iy] = aero.lift([dv_mach,dv_re,dv_aoa,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

    #     ax.plot_surface(MACH, AOA, COEFF, norm=colors.Normalize(vmin=model_range[0],vmax=model_range[1]), cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)

    # plt.xlabel('Mach', fontsize=18)
    # plt.ylabel('AoA', fontsize=18)

    # fig.savefig(os.path.join(project_folder,'fig_2.png'))
    
    # # plt.draw()
    # # plt.show()

    # plt.close(fig)



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
