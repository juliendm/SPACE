#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, glob, re
# from .. import io   as spaceio
# from .  import func as spacefunc
# from ..io import redirect_folder, save_data

import subprocess
import numpy

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

SPACE_RUN = os.environ['SPACE_RUN']

# from optparse import OptionParser
sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE import io   as spaceio
from SPACE.eval import func as spacefunc
from SPACE.eval import design as spacedesign
from SPACE.io import redirect_folder, save_data

from SPACE.util import LHC_unif, DesignVariables
from SPACE.surfpack import Surfpack

# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------

def create_sol():

    desvar = DesignVariables()

    struct_lift_model = Surfpack('STRUCT_LIFT', desvar.ndim-2)
    struct_lift_model.load_model('struct_lift.sps')
    struct_drag_model = Surfpack('STRUCT_DRAG', desvar.ndim-2)
    struct_drag_model.load_model('struct_drag.sps')

    lift_model = Surfpack('LIFT', desvar.ndim)
    lift_model.load_model('../FLIGHT_MECHANICS/model_lift_sup.sps')
    drag_model = Surfpack('DRAG', desvar.ndim)
    drag_model.load_model('../FLIGHT_MECHANICS/model_drag_sup.sps')

    mach_vec = [1.2, 1.5, 1.7, 2.0, 3.0, 5.0, 7.0]
    aoa_range  = [
                      [0.0,15.0],
                      [0.0,17.0],
                      [2.5,19.0],
                      [2.5,20.0],
                      [5.0,25.0],
                      [7.5,35.0],
                      [15.0,40.0]
    ]

    nx = 100

    fig_1 = plt.figure(1)
    fig_2 = plt.figure(2)

    for mach_index in range(len(mach_vec)):

        mach = mach_vec[mach_index]

        aoa_vec = numpy.linspace(aoa_range[mach_index][0], aoa_range[mach_index][1], nx)

        lift_vec = numpy.zeros(len(aoa_vec))
        drag_vec = numpy.zeros(len(aoa_vec))

        struct_lift_vec = numpy.zeros(len(aoa_vec))
        struct_drag_vec = numpy.zeros(len(aoa_vec))

        for aoa_index in range(len(aoa_vec)):

            dvs = desvar.reverse_variable_change([mach,1E7,aoa_vec[aoa_index], 0.0,0.0, -0.5,0.5,0.0,0.0, 0.0,0.0])
            dvs[desvar.rey_index] = 0.0 # Force Reynolds to Reference Trajectory

            dvs_filtered = dvs[0:3] + dvs[5:11]

            lift_vec[aoa_index] = lift_model.eval(dvs)
            drag_vec[aoa_index] = drag_model.eval(dvs)

            struct_lift_vec[aoa_index] = struct_lift_model.eval(dvs_filtered)
            struct_drag_vec[aoa_index] = struct_drag_model.eval(dvs_filtered)

        plt.figure(1)
        plt.plot(aoa_vec,lift_vec,'b')
        plt.plot(aoa_vec,struct_lift_vec,'r')

        plt.figure(2)
        plt.plot(aoa_vec,drag_vec,'b')
        plt.plot(aoa_vec,struct_drag_vec,'r')

    plt.figure(1)
    plt.title('Lift Coefficient', fontsize=18)
    plt.xlabel('AoA [Deg]', fontsize=18)
    plt.legend(mach_vec, title='Mach')

    plt.figure(2)
    plt.title('Drag Coefficient', fontsize=18)
    plt.xlabel('AoA [Deg]', fontsize=18)
    plt.legend(mach_vec, title='Mach')

    fig_1.savefig('fig_lift.png')
    fig_2.savefig('fig_drag.png')

    plt.close(fig_1)
    plt.close(fig_2)




# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':

    create_sol()



