#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
import scipy as sp
from optparse import OptionParser
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

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

    # Mission

    # mission_data = np.loadtxt('output_spaceplane.dat',skiprows=1)
    # mission_phase = mission_data[:,1]
    # mission_altitude = mission_data[:,2]
    # indexes_phase_asc = [n for n,i in enumerate(mission_phase) if i in [1,2]]
    # indexes_phase_dsc = [n for n,i in enumerate(mission_phase) if i in [3,4,5,6]]
    # indexes_phase_noaero = [n for n,i in enumerate(mission_altitude) if i > 80.0]
    # indexes_phase_dsc_noaero = []
    # indexes_phase_dsc_aero = []
    # for index in indexes_phase_dsc:
    #     if index in indexes_phase_noaero:
    #         indexes_phase_dsc_noaero.append(index)
    #     else:
    #         indexes_phase_dsc_aero.append(index)
    # mission_aoa = mission_data[:,9]
    # mission_mach = mission_data[:,26]
    # mission_reynolds = mission_data[:,27]

    # Plot

    nx = 100

    mach_range = [ 0.0, 9.0]
    aoa_range  = [ 0.0,55.0]

    cogs = [0.60*17, 0.63*17]
    cog_z = -0.3

    ndim = 7
    aero_desvar = DesignVariables()

    model_trim_fwd = Surfpack('TRIM_FWD',ndim)
    model_trim_fwd.load_model(os.path.join(project_folder,'model_TRIM_FWD.sps'))
    model_trim_aft = Surfpack('TRIM_AFT',ndim)
    model_trim_aft.load_model(os.path.join(project_folder,'model_TRIM_AFT.sps'))

    model_trimmed_static_margin_aft = Surfpack('TRIMMED_STATIC_MARGIN_AFT',ndim)
    model_trimmed_static_margin_aft.load_model(os.path.join(project_folder,'model_TRIMMED_STATIC_MARGIN_AFT.sps'))
    model_trimmed_static_margin_fwd = Surfpack('TRIMMED_STATIC_MARGIN_FWD',ndim)
    model_trimmed_static_margin_fwd.load_model(os.path.join(project_folder,'model_TRIMMED_STATIC_MARGIN_FWD.sps'))

    for flag in ['TRIM']: # ['TRIMMED_STATIC_MARGIN']: # 

        for shift in [0.0]:


            # # Small
            # dv1 = 0.5 + shift    # Leading Edge Location
            # dv2 = -0.5           # Elevon Length
            # dv3 = -0.5           # Span
            # dv4 = -0.2 + shift   # Hinge Location
            # dv5 = 0.5
            # dv6 = 0.5


            # # Elevon
            dv1 = 0.5 + shift    # Leading Edge Location
            dv2 = 0.5           # Elevon Length
            dv3 = -0.5           # Span
            dv4 = 0.5 + shift   # Hinge Location
            dv5 = 0.5
            dv6 = 0.5

            fig = plt.figure()

            plt.title(flag.replace('_', ' '), fontsize=18)
            plt.xlabel('Mach', fontsize=18)
            plt.ylabel('AoA [Deg]', fontsize=18)

            plt.xlim([0,9])
            plt.ylim([0,55])
            plt.fill_between([0,1,8], [15,15,55], [55,55,55], color='grey', alpha='0.3', lw=0, zorder=3)
            plt.fill_between([1,8,8,9], [0,0,0,0], [0,15,55,55], color='grey', alpha='0.3', lw=0, zorder=3)

            for index,cog_x in enumerate(cogs):

                MACH = np.linspace(mach_range[0], mach_range[1], nx)
                AOA = np.linspace(aoa_range[0], aoa_range[1], nx)
                MACH, AOA = np.meshgrid(MACH, AOA)
                COEFF = MACH*0.0

                for ix in range(nx):
                    for iy in range(nx):


                        aero_dvs = aero_desvar.reverse_variable_change([MACH[ix,iy],1E7,AOA[ix,iy], 0.0,0.0, dv1,dv2,dv3,dv4, dv5,dv6])  ##### MUST DO BETTER for Reynolds... Must be the one from trajectory!!! FOR REYNOLDS
                        aero_dvs[aero_desvar.rey_index] = 0.0 # Force Reynolds to Reference Trajectory
                        perfo_dvs = aero_to_perfo(aero_dvs)

                        if (flag == 'TRIM'):

                            if index == 0:
                                val = model_trim_fwd.eval(perfo_dvs)
                            else:
                                val = model_trim_aft.eval(perfo_dvs)

                            if val < aero_desvar.bf_bound[0]:
                                val = val - aero_desvar.bf_bound[0]
                            elif val > aero_desvar.bf_bound[1]:
                                val = val - aero_desvar.bf_bound[1]
                            else:
                                val = 0.0

                            COEFF[ix,iy] = val*180.0/np.pi

                        elif (flag == 'TRIMMED_STATIC_MARGIN'):

                            if index == 0:
                                COEFF[ix,iy] = model_trimmed_static_margin_fwd.eval(perfo_dvs)
                            else:
                                COEFF[ix,iy] = model_trimmed_static_margin_aft.eval(perfo_dvs)

                        elif (flag == 'K_ALPHA'):

                            COEFF[ix,iy] = 0.0

                        elif (flag == 'EFFICIENCY'):

                            COEFF[ix,iy] = 0.0


                if (flag == 'TRIM'):

                    val = -20.0*0.2
                    levels = [-float("inf"),val];
                    if index == 0:
                        hatches = ['\\\\']
                    else:
                        hatches = ['//']
                    plt.contourf(MACH, AOA, COEFF, levels, colors= 'none', edgecolor='r', hatches=hatches)
                    levels = [val];
                    CS = plt.contour(MACH, AOA, COEFF, levels, colors='k')
                    plt.clabel(CS, inline=1, fontsize=10)

                    val = 30.0*0.2
                    levels = [val,float("inf")];
                    if index == 0:
                        hatches = ['\\\\']
                    else:
                        hatches = ['//']
                    plt.contourf(MACH, AOA, COEFF, levels, colors= 'none', edgecolor='r', hatches=hatches)
                    levels = [val];
                    CS = plt.contour(MACH, AOA, COEFF, levels, colors='k')
                    plt.clabel(CS, inline=1, fontsize=10)

                elif (flag == 'TRIMMED_STATIC_MARGIN'):

                    val = -0.04
                    levels = [-10.0,val];
                    if index == 0:
                        hatches = ['\\\\']
                    else:
                        hatches = ['//']
                    plt.contourf(MACH, AOA, COEFF, levels, colors= 'none', edgecolor='r', hatches=hatches, extend='lower')
                    levels = [val];
                    CS = plt.contour(MACH, AOA, COEFF, levels, colors='k')
                    plt.clabel(CS, inline=1, fontsize=10)

                elif (flag == 'K_ALPHA'):

                    val = 4.0
                    levels = [val,10000.0];
                    if index == 0:
                        hatches = ['\\\\']
                    else:
                        hatches = ['//']
                    plt.contourf(MACH, AOA, COEFF, levels, colors= 'none', edgecolor='r', hatches=hatches, extend='lower')
                    levels = [val];
                    CS = plt.contour(MACH, AOA, COEFF, levels, colors='k')
                    plt.clabel(CS, inline=1, fontsize=10)

                elif (flag == 'EFFICIENCY'):

                    levels = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
                    plt.contour(MACH, AOA, COEFF, levels)


            # plt.plot(mission_mach[indexes_phase_asc], mission_aoa[indexes_phase_asc], color='red')
            # plt.plot(mission_mach[indexes_phase_dsc_noaero], mission_aoa[indexes_phase_dsc_noaero], color='green')
            # plt.plot(mission_mach[indexes_phase_dsc_aero], mission_aoa[indexes_phase_dsc_aero], color='blue')

            fig.savefig(os.path.join(project_folder,'fig_' + flag + '_' + str(shift) + '.png'))
            plt.close(fig)

def aero_to_perfo(dvs):

    return [dvs[0],dvs[1],dvs[2],dvs[5],dvs[6],dvs[7],dvs[8]]


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
