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


from optparse import OptionParser

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

def process_design(regime = 'ON'):

    original = False

    range_number = 300

    if original:
       build_points_folder = 'BUILD_POINTS_ORIGINAL'
       designs_folder = 'DESIGNS_BATCH_1'
    else:
        build_points_folder = 'BUILD_POINTS'
        designs_folder = 'DESIGNS'

    desvar = DesignVariables()

    if regime == 'ON':

        ndim_struct = desvar.ndim_struct_on
        pack_structure = desvar.pack_structure_on

    elif regime == 'OFF':

        ndim_struct = desvar.ndim_struct_off
        pack_structure = desvar.pack_structure_off

    elif regime == 'BOTH':

        ndim_struct = desvar.ndim_struct
        pack_structure = desvar.pack_structure

    threshold = 1.0001




    mass_models = []
    for index in range(132):
        mass_model = Surfpack('MASS_%05d' % (index+1), ndim_struct)
        # if original:
        #     mass_model.load_data(os.path.join(build_points_folder,'build_points_mass_%05d.dat' % (index+1)))
        # else:
        #     mass_model.load_data(os.path.join('BUILD_POINTS_ORIGINAL','enriched_points_mass_%05d.dat' % (index+1))) 
        mass_models.append(mass_model)

    thickness_models = []
    for index in range(132):
        thickness_model = Surfpack('THICKNESS_%05d' % (index+1), ndim_struct)
        # if original:
        #     thickness_model.load_data(os.path.join(build_points_folder,'build_points_thickness_%05d.dat' % (index+1)))
        # else:
        #     thickness_model.load_data(os.path.join('BUILD_POINTS_ORIGINAL','enriched_points_thickness_%05d.dat' % (index+1))) 
        thickness_models.append(thickness_model)

    area_models = []
    for index in range(132):
        area_model = Surfpack('AREA_%05d' % (index+1), ndim_struct)
        # if original:
        #     area_model.load_data(os.path.join(build_points_folder,'build_points_area_%05d.dat' % (index+1)))
        # else:
        #     area_model.load_data(os.path.join('BUILD_POINTS_ORIGINAL','enriched_points_area_%05d.dat' % (index+1))) 
        area_models.append(area_model)




    structure_mass_model = Surfpack('STRUCTURE_MASS', ndim_struct)
    structure_area_model = Surfpack('STRUCTURE_AREA', ndim_struct)
    # if original:
    #     structure_mass_model.load_data(os.path.join(build_points_folder,'build_points_STRUCTURE_MASS.dat'))
    #     structure_area_model.load_data(os.path.join(build_points_folder,'build_points_STRUCTURE_AREA.dat'))
    # else:
    #     structure_mass_model.load_data(os.path.join('BUILD_POINTS_ORIGINAL','enriched_points_STRUCTURE_MASS.dat'))
    #     structure_area_model.load_data(os.path.join('BUILD_POINTS_ORIGINAL','enriched_points_STRUCTURE_AREA.dat'))

    for dsn_index in range(range_number):

        design_folder = os.path.join(designs_folder,'DSN_%03d' % (dsn_index+1))
        history_file_name = os.path.join(design_folder,'STRUCTURE/history_structure.dat')
        mass_file_name = os.path.join(design_folder,'STRUCTURE/lc0_mass_member.dat')
        thickness_file_name = os.path.join(design_folder,'STRUCTURE/lc0_x_final.dat')
        postpro_file_name = os.path.join(design_folder,'STRUCTURE/postpro_load_1.dat')

        if os.path.exists(history_file_name) and os.path.exists(mass_file_name) and os.path.exists(postpro_file_name) and os.path.exists(thickness_file_name):
            history = numpy.loadtxt(history_file_name,delimiter=',',skiprows=1)
            structure_mass = history[-1,1]
            check = history[-1,2:5]
            if (check[0] < threshold) and (check[1] < threshold) and (check[2] < threshold):
                print dsn_index+1,': Success'
                config = spaceio.Config(os.path.join(design_folder,'config_DSN.cfg'))
                dvs = pack_structure(config)
                structure_mass_model.add(dvs, structure_mass)

                mass_data = numpy.loadtxt(mass_file_name)
                for index,mass_model in enumerate(mass_models):
                    mass_model.add(dvs, mass_data[index])

                thickness_data = numpy.loadtxt(thickness_file_name)
                for index,thickness_model in enumerate(thickness_models):
                    thickness_model.add(dvs, thickness_data[index])

                structure_area = 0
                for index,area_model in enumerate(area_models):
                    area = mass_data[index]/thickness_data[index]/float(config.MATERIAL_DENSITY)
                    area_model.add(dvs, area)
                    structure_area += area

                structure_area_model.add(dvs, structure_area)


            else:
                print dsn_index+1, ': Warning: ',check[0],check[1],check[2]
        else:
            print 'Missing:', dsn_index+1

    structure_mass_model.save_data(os.path.join(build_points_folder, 'enriched_points_STRUCTURE_MASS.dat'))
    structure_area_model.save_data(os.path.join(build_points_folder, 'enriched_points_STRUCTURE_AREA.dat'))

    for index,mass_model in enumerate(mass_models):
        print 'Save MASS: ', index
        mass_model.save_data(os.path.join(build_points_folder, 'enriched_points_mass_%05d.dat' % (index+1)))

    for index,thickness_model in enumerate(thickness_models):
        print 'Save THICKNESS: ', index
        thickness_model.save_data(os.path.join(build_points_folder, 'enriched_points_thickness_%05d.dat' % (index+1)))

    for index,area_model in enumerate(area_models):
        print 'Save AREA: ', index
        area_model.save_data(os.path.join(build_points_folder, 'enriched_points_area_%05d.dat' % (index+1)))

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-r", "--regime", dest="regime", default="ON",
                      help="regime", metavar="REGIME")


    (options, args)=parser.parse_args()

    process_design( options.regime )


