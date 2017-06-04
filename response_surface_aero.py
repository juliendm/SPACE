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
        config = SPACE.io.Config(filename)
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    print '%d design(s) so far' % len(project.designs)

    konfig = copy.deepcopy(config)
    konfig.NUMBER_PART = partitions

    # Design Variables

    desvar = DesignVariables()

    if initiate:

        nd = 10*desvar.ndim

        X = LHC_unif(desvar.XB_SUP,nd)

        lift_model = Surfpack('LIFT',desvar.ndim)
        drag_model = Surfpack('DRAG',desvar.ndim)
        moment_y_model = Surfpack('MOMENT_Y',desvar.ndim)

        for index in range(0,len(X)):

            dvs = X[index]
            desvar.unpack(konfig, dvs)

            lift_model.add(dvs,project.func('LIFT',konfig))
            drag_model.add(dvs,project.func('DRAG',konfig))
            moment_y_model.add(dvs,project.func('MOMENT_Y',konfig))

        lift_model.save_data(os.path.join(project_folder,'build_points_lift.dat'))
        drag_model.save_data(os.path.join(project_folder,'build_points_drag.dat'))
        moment_y_model.save_data(os.path.join(project_folder,'build_points_moment_y.dat'))

    else:

        na = 500

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

        for ite in range(0, na): # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if ite%3 == 0:
                dvs = lift_model.max_variance(desvar.XB_SUP)
            elif ite%3 == 1:
                dvs = drag_model.max_variance(desvar.XB_SUP)
            else:
                dvs = moment_y_model.max_variance(desvar.XB_SUP)

            desvar.unpack(konfig, dvs)

            print '-------------------------------'
            print dvs
            print '-------------------------------'

            proc = project.func('AERODYNAMICS', konfig)
            new_design_container = project.designs[-1]

            proc.wait()
            project.compile_designs()

            new_design = new_design_container.design
            if not new_design is None:
                new_dvs = desvar.pack(new_design.config)
                new_funcs = new_design.funcs
                if hasattr(new_funcs,'LIFT') and hasattr(new_funcs,'DRAG') and hasattr(new_funcs,'MOMENT_Y'):
                    lift_model.add(new_dvs,new_funcs.LIFT)
                    drag_model.add(new_dvs,new_funcs.DRAG)
                    moment_y_model.add(new_dvs,new_funcs.MOMENT_Y)

            # Build Model
            lift_model.build('kriging')
            drag_model.build('kriging')
            moment_y_model.build('kriging')

            lift_model.save_data(os.path.join(project_folder,'enriched_points_lift.dat'))
            drag_model.save_data(os.path.join(project_folder,'enriched_points_drag.dat'))
            moment_y_model.save_data(os.path.join(project_folder,'enriched_points_moment_y.dat'))

        lift_model.save_model(os.path.join(project_folder,'model_lift.sps'))
        drag_model.save_model(os.path.join(project_folder,'model_drag.sps'))
        moment_y_model.save_model(os.path.join(project_folder,'model_moment_y.sps'))

#: response_surface()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
