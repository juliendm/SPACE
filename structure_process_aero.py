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

def process_design():

    range_number = 300

    desvar = DesignVariables()

    build_points_folder = 'BUILD_POINTS'

    lift_model = Surfpack('LIFT', desvar.ndim-2)
    lift_model.load_data(os.path.join(build_points_folder,'build_points_struct_lift_model.dat'))
    drag_model = Surfpack('DRAG', desvar.ndim-2)
    drag_model.load_data(os.path.join(build_points_folder,'build_points_struct_drag_model.dat'))


    threshold = 1.0001

    for dsn_index in range(range_number):

        design_folder = 'DESIGNS/DSN_%03d' % (dsn_index+1)

        print design_folder

        postpro_file_name = os.path.join(design_folder,'STRUCTURE/postpro_load_1.dat')


        history_file_name = os.path.join(design_folder,'STRUCTURE/history_structure.dat')
        mass_file_name = os.path.join(design_folder,'STRUCTURE/lc0_mass_member.dat')
        postpro_file_name = os.path.join(design_folder,'STRUCTURE/postpro_load_1.dat')

        if os.path.exists(history_file_name) and os.path.exists(mass_file_name) and os.path.exists(postpro_file_name):

            history = numpy.loadtxt(history_file_name,delimiter=',',skiprows=1)
            structure_mass = history[-1,1]
            check = history[-1,2:5]
            if (check[0] < threshold) and (check[1] < threshold) and (check[2] < threshold):

                config = spaceio.Config(os.path.join(design_folder,'config_DSN.cfg'))
                dvs = desvar.pack(config)
                dvs_filtered = dvs[0:3] + dvs[5:11]

                print dvs_filtered

                with open(postpro_file_name) as fp:
                    for i, line in enumerate(fp):
                        if i == 13:
                            data = line.split(':')[-1].split(',')
                            lift = float(data[0])
                            drag = float(data[1])
                            print lift, drag
                        elif i > 13:
                            break

                lift_model.add(dvs_filtered, lift)
                drag_model.add(dvs_filtered, drag)


    lift_model.save_data(os.path.join(build_points_folder,'enriched_points_struct_lift_model.dat'))
    drag_model.save_data(os.path.join(build_points_folder,'enriched_points_struct_drag_model.dat'))

    lift_model.build('kriging')
    drag_model.build('kriging')

    lift_model.save_model('struct_lift.sps')
    drag_model.save_model('struct_drag.sps')


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':

    process_design()



