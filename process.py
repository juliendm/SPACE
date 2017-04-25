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
                      
    (options, args)=parser.parse_args()

    process( options.project_folder )

#: main()

def process( project_folder ):

    project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
    project.compile_designs()

    desvar = DesignVariables()

    lift_model = Surfpack('LIFT',desvar.ndim)
    drag_model = Surfpack('DRAG',desvar.ndim)
    moment_y_model = Surfpack('MOMENT_Y',desvar.ndim)

    for design_container in project.designs:
        design = design_container.design
        if not design is None:
            dvs = desvar.pack(design.config)
            funcs = design.funcs
            if hasattr(funcs,'LIFT') and hasattr(funcs,'DRAG') and hasattr(funcs,'MOMENT_Y'):
                lift_model.add(dvs,funcs.LIFT)
                drag_model.add(dvs,funcs.DRAG)
                moment_y_model.add(dvs,funcs.MOMENT_Y)
                print 'done ' + design.folder
            else:
                print 'missing ' + design.folder
        else:
            print 'missing ' + design.folder

    lift_model.save_data(os.path.join(project_folder,'build_points_lift.dat'))
    drag_model.save_data(os.path.join(project_folder,'build_points_drag.dat'))
    moment_y_model.save_data(os.path.join(project_folder,'build_points_moment_y.dat'))

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()

