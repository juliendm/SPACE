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
    
    aero = spacemodel.AeroModel(konfig,project_folder)

    # print aero.static_margin([8.000000e+00,   1.000000e+00,  -1.000000e+00,   2.059314e-01,   5.000000e-01,   5.000000e-01,  -5.000000e-01,   5.000000e-01,  -2.000000e-01,  -5.000000e-01,  -2.000000e-01],11.0)
    # print aero.moment_y([5.710816e+00,   1.000000e+00,  -1.000000e+00,  -2.500000e-01,  -3.500000e-01,  -5.000000e-01,  -5.000000e-01,  -5.000000e-01,  -2.000000e-01,  -5.000000e-01,  -2.000000e-01],11.0)
    # print aero.moment_y([3.665795e+00,   1.000000e+00,  -1.000000e+00,   5.000000e-01,   5.000000e-01,  -5.000000e-01,  -5.000000e-01,   5.000000e-01,  -2.000000e-01,  -5.000000e-01,  -2.000000e-01],11.0)

    # print aero.max_trimmed_efficiency([0.5,0.0,0.0,  0.0,0.0,  -0.5,0.5,0.5,0.5,  0.5,0.5],10.0)

    # print aero.moment_y([5.623184e+00,   1.000000e+00,   1.000000e+00,   9.860460e-02,  -2.500000e-01,   5.000000e-01,  -5.000000e-01,   5.000000e-01,   5.000000e-01,   5.000000e-01,  -2.000000e-01], 11.0)
    
    # dvs = [2.377544e+00,  -5.045048e-01,   8.128633e-01,   2.803823e-01 ,  3.110862e-01,  -4.362524e-01 ,  1.069665e-01 , -1.411665e-01 , -7.990462e-02 ,  3.137270e-01,   1.857822e-01]
    # print aero.drag(dvs) #   3.007715e-01
    # print aero.lift(dvs) #   3.007715e-01
    # print aero.moment_y(dvs,11.0) #   3.007715e-01



    # Plot

    surface = True

    if surface:

        flag = 'TRIM'

        nx = 100

        mach_range = [[ 0.0, 8.0]]
        aoa_range  = [[ 0.0,55.0]]

        cogs = [10.5] # [10.1,10.3,10.5,10.7,10.9] # 9.6,9.7,9.8,9.9,
        cog_z = -0.5

        for cog_x in cogs:

            fig = plt.figure()

            for index in range(len(mach_range)):

                MACH = np.linspace(mach_range[index][0], mach_range[index][1], nx)
                AOA = np.linspace(aoa_range[index][0], aoa_range[index][1], nx)
                MACH, AOA = np.meshgrid(MACH, AOA)
                COEFF = MACH*0.0

                for ix in range(nx):
                    for iy in range(nx):

                        print ix, iy

                        dvs = aero.desvar.reverse_variable_change([MACH[ix,iy],1E7,AOA[ix,iy], 0.0,0.0, -0.5,0.5,0.5,0.5, 0.5,0.5])  ##### MUST DO BETTER for Reynolds... Must be the one from trajectory!!! FOR REYNOLDS
                        dvs[aero.desvar.rey_index] = 0.0 # Force Reynolds to Reference Trajectory

                        if (dvs[aero.desvar.aoa_index] <= 1.0 and dvs[aero.desvar.aoa_index] >= -1.0):

                            if (flag == 'DRAG'):

                                COEFF[ix,iy] = aero.drag(dvs)

                            elif (flag == 'PITCH'):

                                COEFF[ix,iy] = aero.moment_y(dvs,cog_x,cog_z)

                            elif (flag == 'TRIM'):

                                COEFF[ix,iy] = aero.trim(dvs,cog_x,cog_z)[aero.desvar.el_index]*180.0/np.pi
                            
                            elif (flag == 'STATIC_MARGIN'):

                                trimmed_dvs = aero.trim(dvs,cog_x,cog_z)
                                if (trimmed_dvs[aero.desvar.el_index]*180.0/np.pi < 30.0*0.3) and (trimmed_dvs[aero.desvar.el_index]*180.0/np.pi > -20.0*0.3):
                                    COEFF[ix,iy] = aero.static_margin(trimmed_dvs,cog_x,cog_z)
                                else:
                                    COEFF[ix,iy] = -float("inf") #float('nan')

                            elif (flag == 'K_ALPHA'):

                                trimmed_dvs = aero.trim(dvs,cog_x,cog_z)
                                #if (trimmed_dvs[aero.desvar.el_index]*180.0/np.pi < 30.0*0.3) and (trimmed_dvs[aero.desvar.el_index]*180.0/np.pi > -20.0*0.3):
                                COEFF[ix,iy] = aero.k_alpha(trimmed_dvs,cog_x,cog_z)
                                #else:
                                #    COEFF[ix,iy] = -float("inf") #float('nan')
                            
                            elif (flag == 'EFFICIENCY'):

                                COEFF[ix,iy] = aero.max_trimmed_efficiency(dvs,cog_x,cog_z)

                        else:

                            COEFF[ix,iy] = float('NaN')

                if (flag == 'DRAG'):

                    CS = plt.contour(MACH, AOA, COEFF, [0.0, 0.1, 0.2, 0.3])

                elif (flag == 'PITCH'):

                    CS = plt.contour(MACH, AOA, COEFF, [-0.05,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.05])

                elif (flag == 'TRIM'):

                    levels = [-20.0*0.3, -20.0*0.2, -20.0*0.1, 0.0, 30.0*0.1, 30.0*0.2, 30.0*0.3]; colors = ('r', 'g', 'r')
                    CS = plt.contour(MACH, AOA, COEFF, levels)

                elif (flag == 'STATIC_MARGIN'):

                    levels = [-0.05, -0.04, -0.03]; colors = ('r', 'g')
                    CS = plt.contour(MACH, AOA, COEFF, levels)

                elif (flag == 'K_ALPHA'):

                    levels = [1.0, 2.0, 3.0, 4.0]
                    CS = plt.contour(MACH, AOA, COEFF, levels)

                elif (flag == 'EFFICIENCY'):

                    levels = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
                    CS = plt.contour(MACH, AOA, COEFF, levels)

                plt.clabel(CS, inline=1, fontsize=10)

            plt.xlabel('Mach', fontsize=18)
            plt.ylabel('AoA', fontsize=18)

            plt.ylim([0,55])
            plt.fill_between([0,1], [15,15], [55,55], color='grey', alpha='0.5', zorder=3)
            plt.fill_between([1,8], [15,55], [55,55], color='grey', alpha='0.5', zorder=3)
            plt.fill_between([1,8], [0,0], [0,15], color='grey', alpha='0.5', zorder=3)

            fig.savefig(os.path.join(project_folder,'fig_' + str(cog_x) + '.png'))
            plt.close(fig)


    else:

        mach_vec = [0.3, 0.7, 0.9, 1.2, 1.5, 1.7, 2.0, 3.0, 5.0, 7.0]
        aoa_range  = [[0.0,15.0],
                      [0.0,15.0],
                      [0.0,15.0],
                      [0.0,15.0],
                      [0.0,17.0],
                      [2.5,19.0],
                      [2.5,20.0],
                      [5.0,25.0],
                      [7.5,35.0],
                      [15.0,40.0]]

        cog = 9.7

        nx = 100

        fig_1 = plt.figure(1)
        fig_2 = plt.figure(2)
        fig_3 = plt.figure(3)

        for mach_index in range(len(mach_vec)):

            mach = mach_vec[mach_index]

            aoa_vec = np.linspace(aoa_range[mach_index][0], aoa_range[mach_index][1], nx)
            lift_vec = np.zeros(len(aoa_vec))
            drag_vec = np.zeros(len(aoa_vec))
            moment_y_vec = np.zeros(len(aoa_vec))

            for aoa_index in range(len(aoa_vec)):

                dvs = aero.desvar.reverse_variable_change([mach,1E7,aoa_vec[aoa_index], 0.0,0.0, -0.5,0.5,0.0,0.0, 0.0,0.0])
                dvs[aero.desvar.rey_index] = 0.0 # Force Reynolds to Reference Trajectory

                lift_vec[aoa_index] = aero.lift(dvs)
                drag_vec[aoa_index] = aero.drag(dvs)
                moment_y_vec[aoa_index] = aero.moment_y(dvs,9.7,-1.0)

            plt.figure(1)
            plt.plot(aoa_vec,lift_vec)

            plt.figure(2)
            plt.plot(aoa_vec,drag_vec)

            plt.figure(3)
            plt.plot(aoa_vec,moment_y_vec)

        plt.figure(1)
        plt.legend(mach_vec)

        plt.figure(2)
        plt.legend(mach_vec)

        plt.figure(3)
        plt.legend(mach_vec)

        fig_1.savefig(os.path.join(project_folder,'fig_lift.png'))
        fig_2.savefig(os.path.join(project_folder,'fig_drag.png'))
        fig_3.savefig(os.path.join(project_folder,'fig_moment_y.png'))

        plt.close(fig_1)
        plt.close(fig_2)
        plt.close(fig_3)

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
