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

# rice

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

    mission_data = np.loadtxt('output_spaceplane.dat',skiprows=1)
    mission_phase = mission_data[:,1]
    mission_altitude = mission_data[:,2]
    indexes_phase_asc = [n for n,i in enumerate(mission_phase) if i in [1,2]]
    indexes_phase_dsc = [n for n,i in enumerate(mission_phase) if i in [3,4,5,6]]
    indexes_phase_noaero = [n for n,i in enumerate(mission_altitude) if i > 80.0]
    indexes_phase_dsc_noaero = []
    indexes_phase_dsc_aero = []
    for index in indexes_phase_dsc:
        if index in indexes_phase_noaero:
            indexes_phase_dsc_noaero.append(index)
        else:
            indexes_phase_dsc_aero.append(index)
    mission_aoa = mission_data[:,9]
    mission_mach = mission_data[:,26]
    mission_reynolds = mission_data[:,27]

    # Plot

    change_var = False

    surface = False
    case = 'CASE_1'
    compute_COEFF = False

    curves = True

    if change_var:

        fig = plt.figure()
        plt.xlabel('Mach', fontsize=18)
        plt.ylabel('AoA [Deg]', fontsize=18)
        plt.xlim([0,9])
        plt.ylim([0,55])
        plt.fill_between([0,1,8], [15,15,55], [55,55,55], color='grey', alpha='0.3', lw=0, zorder=3)
        plt.fill_between([1,8,8,9], [0,0,0,0], [0,15,55,55], color='grey', alpha='0.3', lw=0, zorder=3)
        fig.savefig(os.path.join(project_folder,'fig_change_var_aoa.png'))
        plt.close(fig)

        fig = plt.figure()
        plt.xlabel('Mach', fontsize=18)
        plt.ylabel('Reynolds', fontsize=18)
        plt.xlim([0,9])
        fill_mach = np.linspace(0,9,100)
        fill_reynolds = 10.0**(-3.0/8.0*fill_mach+7.0+1.0)
        max_y = max(fill_reynolds)
        plt.fill_between(fill_mach, fill_reynolds, np.max(fill_reynolds)*np.ones(len(fill_reynolds)), color='grey', alpha='0.3', zorder=3)
        fill_reynolds = 10.0**(-3.0/8.0*fill_mach+7.0-1.0)
        min_y = min(fill_reynolds)
        plt.fill_between(fill_mach, np.min(fill_reynolds)*np.ones(len(fill_reynolds)), fill_reynolds, color='grey', alpha='0.3', zorder=3)
        plt.ylim([min_y,max_y])
        plt.yscale('log')
        fig.savefig(os.path.join(project_folder,'fig_change_var_reynolds.png'))
        plt.close(fig)


    if surface:

        nx = 50

        mach_range = [ 0.0, 9.0]
        aoa_range  = [ 0.0,55.0]

        cogs = [0.60*17, 0.63*17]
        cog_z = -0.3

        for flag in ['K_ALPHA','TRIM']:

            for shift in [0.0]:

                if case == 'CASE_1':

                    dv1 = 0.5 + shift   # Leading Edge Location
                    dv2 = -0.5          # Elevon Length
                    dv3 = -0.5          # Span
                    dv4 = -0.2 + shift  # Hinge Location

                elif case == 'CASE_3':

                    dv1 = 0.5 + shift   # Leading Edge Location
                    dv2 = 0.48          # Elevon Length
                    dv3 = -0.5          # Span
                    dv4 = 0.39 + shift  # Hinge Location

                else:

                    dv1 = 0.0; dv2 = 0.0; dv3 = 0.0; dv4 = 0.0

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

                    if compute_COEFF:

                        for ix in range(nx):
                            for iy in range(nx):

                                print cog_x, ix, iy

                                dvs = aero.desvar.reverse_variable_change([MACH[ix,iy],1E7,AOA[ix,iy], 0.0,0.0, dv1,dv2,dv3,dv4, dv5,dv6])  ##### MUST DO BETTER for Reynolds... Must be the one from trajectory!!! FOR REYNOLDS
                                dvs[aero.desvar.rey_index] = 0.0 # Force Reynolds to Reference Trajectory

#                               if (dvs[aero.desvar.aoa_index] <= 1.0 and dvs[aero.desvar.aoa_index] >= -1.0):
                                if 1: #(dvs[aero.desvar.aoa_index] <= 3.0 and dvs[aero.desvar.aoa_index] >= -3.0):

                                    if (flag == 'DRAG'):

                                        COEFF[ix,iy] = aero.drag(dvs)

                                    elif (flag == 'PITCH'):

                                        COEFF[ix,iy] = aero.moment_y(dvs,cog_x,cog_z)

                                    elif (flag == 'TRIM'):

                                        COEFF[ix,iy] = aero.trim(dvs,cog_x,cog_z)[aero.desvar.el_index]*180.0/np.pi

                                    elif (flag == 'TRIMMED_STATIC_MARGIN'):

                                        trimmed_dvs = aero.trim(dvs,cog_x,cog_z)
                                        if (trimmed_dvs[aero.desvar.el_index]*180.0/np.pi < 30.0*0.3) and (trimmed_dvs[aero.desvar.el_index]*180.0/np.pi > -20.0*0.3):
                                            COEFF[ix,iy] = aero.static_margin(trimmed_dvs,cog_x,cog_z)
                                        else:
                                            COEFF[ix,iy] = -10.0

                                    elif (flag == 'K_ALPHA'):

                                        trimmed_dvs = aero.trim(dvs,cog_x,cog_z)
                                        if (trimmed_dvs[aero.desvar.el_index]*180.0/np.pi < 30.0*0.3) and (trimmed_dvs[aero.desvar.el_index]*180.0/np.pi > -20.0*0.3):
                                            COEFF[ix,iy] = aero.k_alpha(trimmed_dvs,cog_x,cog_z)
                                        else:
                                            COEFF[ix,iy] = 10000.0

                                    elif (flag == 'EFFICIENCY'):

                                        COEFF[ix,iy] = aero.max_trimmed_efficiency(dvs,cog_x,cog_z)

                                else:

                                    COEFF[ix,iy] = float('NaN')

                        np.save('%s_%d_%s' % (flag,index,case), COEFF)

                    else:

                        COEFF = np.load(os.path.join(project_folder,'%s_%d_%s.npy' % (flag,index,case)))


                    if (flag == 'DRAG'):

                        plt.contour(MACH, AOA, COEFF, [0.0, 0.1, 0.2, 0.3])

                    elif (flag == 'PITCH'):

                        plt.contour(MACH, AOA, COEFF, [-0.05,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.05])

                    elif (flag == 'TRIM'):

                        # levels = [-20.0*0.3, , -20.0*0.1, 0.0, 30.0*0.1, 30.0*0.2, 30.0*0.3]; colors = ('r', 'g', 'r')
                        # plt.contour(MACH, AOA, COEFF, levels)

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

                        # levels = [1.0, 2.0, 3.0, 4.0]
                        # plt.contour(MACH, AOA, COEFF, levels)

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


                plt.plot(mission_mach[indexes_phase_asc], mission_aoa[indexes_phase_asc], color='red')
                plt.plot(mission_mach[indexes_phase_dsc_noaero], mission_aoa[indexes_phase_dsc_noaero], color='green')
                plt.plot(mission_mach[indexes_phase_dsc_aero], mission_aoa[indexes_phase_dsc_aero], color='blue')

                fig.savefig(os.path.join(project_folder,'fig_' + flag + '_' + str(shift) + '.png'))
                plt.close(fig)


        # Reynolds = 10.0**(-3.0/8.0*dv_mach+7.0+dv_rey) #

        fig = plt.figure()

        fill_mach = np.linspace(0,8,100)

        fill_reynolds = 10.0**(-3.0/8.0*fill_mach+7.0+1.0)
        plt.fill_between(fill_mach, fill_reynolds, np.max(fill_reynolds)*np.ones(len(fill_reynolds)), color='grey', alpha='0.5', zorder=3)

        fill_reynolds = 10.0**(-3.0/8.0*fill_mach+7.0-1.0)
        plt.fill_between(fill_mach, np.min(fill_reynolds)*np.ones(len(fill_reynolds)), fill_reynolds, color='grey', alpha='0.5', zorder=3)

        plt.semilogy(mission_mach[indexes_phase_asc], mission_reynolds[indexes_phase_asc], color='red')
        plt.semilogy(mission_mach[indexes_phase_dsc_noaero], mission_reynolds[indexes_phase_dsc_noaero], color='green')
        plt.semilogy(mission_mach[indexes_phase_dsc_aero], mission_reynolds[indexes_phase_dsc_aero], color='blue')

        fig.savefig(os.path.join(project_folder,'fig_reynolds.png'))
        plt.close(fig)

    if curves:

        print 'fig_lift.png', 'fig_drag.png', 'fig_moment_y.png'

        mach_vec = [0.3, 0.7, 0.9, 1.2, 1.5, 1.7, 2.0, 3.0, 5.0, 7.0, 8.0]
        aoa_range  = [[0.0,15.0],
                      [0.0,15.0],
                      [0.0,15.0],
                      [0.0,15.0],
                      [0.0,17.0],
                      [2.5,19.0],
                      [2.5,20.0],
                      [5.0,25.0],
                      [7.5,35.0],
                      [15.0,40.0],[15.0,40.0]]

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
        plt.title('Lift Coefficient', fontsize=18)
        plt.xlabel('AoA [Deg]', fontsize=18)
        plt.legend(mach_vec, title='Mach')

        plt.figure(2)
        plt.title('Drag Coefficient', fontsize=18)
        plt.xlabel('AoA [Deg]', fontsize=18)
        plt.legend(mach_vec, title='Mach')

        plt.figure(3)
        plt.title('Pitch Moment Coefficient', fontsize=18)
        plt.xlabel('AoA [Deg]', fontsize=18)
        plt.legend(mach_vec, title='Mach')

        fig_1.savefig(os.path.join(project_folder,'fig_lift.png'))
        fig_2.savefig(os.path.join(project_folder,'fig_drag.png'))
        fig_3.savefig(os.path.join(project_folder,'fig_moment_y.png'))

        plt.close(fig_1)
        plt.close(fig_2)
        plt.close(fig_3)


        print 'fig_trajectory.png'

        matplotlib.rcParams['figure.figsize'] = 20, 20

        fig = plt.figure()

        ax = plt.subplot(5,2,1)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Altitute [km]', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,2], color='k')

        ax = plt.subplot(5,2,2)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Velocity [m/s]', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,5], color='k')

        ax = plt.subplot(5,2,3)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Mass [kg]', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,8], color='k')

        ax = plt.subplot(5,2,4)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Thrust [N]', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,20], color='k')

        ax = plt.subplot(5,2,5)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Dynamic Pressure [Pa]', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,15], color='k')

        ax = plt.subplot(5,2,6)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Acceleration', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,17], 'k', label='nx')
        ax.plot(mission_data[:,0], mission_data[:,19], 'k--', label='nz')
        ax.legend(fontsize=18)

        ax = plt.subplot(5,2,7)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Heat Flux [W/m^2]', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,16], color='k')

        ax = plt.subplot(5,2,8)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Flight Path Angle [Deg]', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,6], color='k')

        ax = plt.subplot(5,2,9)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('AoA [Deg]', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,9], color='k')

        ax = plt.subplot(5,2,10)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Coefficients', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,21], 'k', label='Drag')
        ax.plot(mission_data[:,0], mission_data[:,22], 'k--', label='Lift')
        ax.legend(fontsize=18)

        fig.savefig(os.path.join(project_folder,'fig_trajectory.png'), bbox_inches='tight')
        plt.close(fig)


        print 'fig_load_cases.png'

        matplotlib.rcParams['figure.figsize'] = 20, 6

        load_case_1 = 20
        load_case_2 = 220
        load_case_3 = 470

        fig = plt.figure()

        ax = plt.subplot(1,2,1)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Dynamic Pressure [Pa]', fontsize=18)
        ax.fill_between([load_case_1-15, load_case_1+15], [0,0], [14000,14000], color='red', alpha='0.3', lw=0, zorder=3)
        ax.fill_between([load_case_2-15, load_case_2+15], [0,0], [14000,14000], color='green', alpha='0.3', lw=0, zorder=3)
        ax.fill_between([load_case_3-15, load_case_3+15], [0,0], [14000,14000], color='blue', alpha='0.3', lw=0, zorder=3)
        ax.plot(mission_data[:,0], mission_data[:,15], color='k')

        ax = plt.subplot(1,2,2)
        ax.set_xlabel('Time [s]', fontsize=18)
        ax.set_ylabel('Acceleration', fontsize=18)
        ax.plot(mission_data[:,0], mission_data[:,17], 'k', label='nx')
        ax.plot(mission_data[:,0], mission_data[:,19], 'k--', label='nz')
        ax.fill_between([load_case_1-15, load_case_1+15], [-5,-5], [3,3], color='red', alpha='0.3', lw=0, zorder=3)
        ax.fill_between([load_case_2-15, load_case_2+15], [-5,-5], [3,3], color='green', alpha='0.3', lw=0, zorder=3)
        ax.fill_between([load_case_3-15, load_case_3+15], [-5,-5], [3,3], color='blue', alpha='0.3', lw=0, zorder=3)
        ax.legend(fontsize=18)

        fig.savefig(os.path.join(project_folder,'fig_load_cases.png'), bbox_inches='tight')
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
