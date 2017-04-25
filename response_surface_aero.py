#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np
import pwd
import pickle

from optparse import OptionParser

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack
from SPACE.util import LHC_unif, DesignVariables

sys.path.append(os.environ['SU2_RUN'])
import SU2

def main():

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--project", dest="project_folder",
                      help="project folder", metavar="PROJECT_FOLDER")
    parser.add_option("-r", "--regime", dest="regime", default="SUP",
                      help="regime", metavar="REGIME")
    parser.add_option("-i", "--initiate", dest="initiate", default="False",
                      help="initiate", metavar="INITIATE")
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")


    (options, args)=parser.parse_args()
    options.initiate = options.initiate.upper() == 'TRUE'
    options.partitions = int( options.partitions )

    response_surface( options.filename       ,
                      options.project_folder ,
                      options.regime         ,
                      options.initiate       ,
                      options.partitions  )

#: main()

def response_surface( filename          ,
                      project_folder    ,
                      regime = 'SUP'    ,
                      initiate = False  ,
                      partitions  = 0  ):

    # Project

    if os.path.exists(project_folder):
        project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
        project.compile_designs()
        config = project.config
    else:
        config = SPACE.io.Config(options.filename)
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    print '%d design(s) so far' % len(project.designs)

    # Design Variables

    desvar = DesignVariables(regime)

    if initiate:

        nd = 10*desvar.ndim

        X = LHC_unif(desvar.XB,nd)

        lift_model = Surfpack('LIFT',desvar.ndim)
        drag_model = Surfpack('DRAG',desvar.ndim)
        moment_y_model = Surfpack('MOMENT_Y',desvar.ndim)

        for index in range(0,len(X)):

            dvs = X[index]

            konfig = copy.deepcopy(config)
            desvar.unpack(konfig, dvs)

            lift_model.add(dvs,project.func('LIFT',konfig))
            drag_model.add(dvs,project.func('DRAG',konfig))
            moment_y_model.add(dvs,project.func('MOMENT_Y',konfig))

        lift_model.save_data(os.path.join(project_folder,'build_points_lift.dat'))
        drag_model.save_data(os.path.join(project_folder,'build_points_drag.dat'))
        moment_y_model.save_data(os.path.join(project_folder,'build_points_moment_y.dat'))

    else:

        na = 30

        lift_model = Surfpack('LIFT',desvar.ndim)
        lift_model.load_data(os.path.join(project_folder,'build_points_lift.dat'))

        drag_model = Surfpack('DRAG',desvar.ndim)
        drag_model.load_data(os.path.join(project_folder,'build_points_drag.dat'))

        moment_y_model = Surfpack('MOMENT_Y',desvar.ndim)
        moment_y_model.load_data(os.path.join(project_folder,'build_points_moment_y.dat'))

        # Build Model
        lift_model.build('kriging')
        drag_model.build('kriging')
        moment_y_model.build('kriging')

        for ite in range(na):

            dvs = lift_model.max_variance(desvar.XB)

            print '-------------------------------'
            print dvs
            print '-------------------------------'

            break

            # konfig = copy.deepcopy(config)
            # apply_dvs(konfig, dvs)

            # lift_model.add(dvs,project.func('LIFT',konfig))
            # drag_model.add(dvs,project.func('DRAG',konfig))
            # moment_y_model.add(dvs,project.func('MOMENT_Y',konfig))

            # # Build Model
            # lift_model.build('kriging')
            # drag_model.build('kriging')
            # moment_y_model.build('kriging')

            # lift_model.save_data(os.path.join(project_folder,'enriched_points_lift.dat'))
            # drag_model.save_data(os.path.join(project_folder,'enriched_points_drag.dat'))
            # moment_y_model.save_data(os.path.join(project_folder,'enriched_points_moment_y.dat'))

        # lift_model.save_model(os.path.join(project_folder,'lift.sps'))
        # drag_model.save_model(os.path.join(project_folder,'drag.sps'))
        # moment_y_model.save_model(os.path.join(project_folder,'moment_y.sps'))

#: response_surface()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
