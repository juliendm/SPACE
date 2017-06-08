#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
from optparse import OptionParser
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack
from SPACE.util import DesignVariables

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-p", "--project", dest="project_folder",
                      help="project folder", metavar="PROJECT_FOLDER")
    parser.add_option("-r", "--regime", dest="regime", default="SUP",
                      help="regime", metavar="REGIME")
                      
    (options, args)=parser.parse_args()

    process_project( options.project_folder ,
             options.regime )

#: main()

def process_project( project_folder   ,
                     regime = 'SUP' ):

    project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
    project.compile_designs()

    desvar = DesignVariables(regime)

    lift_model = Surfpack('LIFT',desvar.ndim)
    drag_model = Surfpack('DRAG',desvar.ndim)
    force_z_model = Surfpack('FORCE_Z',desvar.ndim)
    moment_y_model = Surfpack('MOMENT_Y',desvar.ndim)

    for design_container in project.designs:
        design = design_container.design
        if not design is None:
            dvs = desvar.pack(design.config)
            funcs = design.funcs
            if hasattr(funcs,'LIFT') and hasattr(funcs,'DRAG') and hasattr(funcs,'MOMENT_Y'):
                lift_model.add(dvs,funcs.LIFT)
                drag_model.add(dvs,funcs.DRAG)
                force_z_model.add(dvs,funcs.FORCE_Z)
                moment_y_model.add(dvs,funcs.MOMENT_Y)
#                print 'done ' + design.folder
            else:
                print 'missing ' + design.folder
        else:
            print 'missing design'

    # Save Data
    lift_model.save_data(os.path.join(project_folder,'build_points_lift.dat'))
    drag_model.save_data(os.path.join(project_folder,'build_points_drag.dat'))
    force_z_model.save_data(os.path.join(project_folder,'build_points_force_z.dat'))
    moment_y_model.save_data(os.path.join(project_folder,'build_points_moment_y.dat'))

    # Build Model
    lift_model.build('kriging')
    drag_model.build('kriging')
    force_z_model.build('kriging')
    moment_y_model.build('kriging')

    # Save Model
    lift_model.save_model(os.path.join(project_folder,'model_lift.sps'))
    drag_model.save_model(os.path.join(project_folder,'model_drag.sps'))
    force_z_model.save_model(os.path.join(project_folder,'model_force_z.sps'))
    moment_y_model.save_model(os.path.join(project_folder,'model_moment_y.sps'))

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()

