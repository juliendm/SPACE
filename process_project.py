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
    parser.add_option("-m", "--max_dsn", dest="max_dsn", default="ALL",
                      help="max design", metavar="MAX_DESIGN")
                      
    (options, args)=parser.parse_args()

    process_project( options.project_folder, options.max_dsn )

#: main()

def process_project( project_folder, max_dsn ):

    desvar = DesignVariables()

    lift_model = Surfpack('LIFT',desvar.ndim)
    drag_model = Surfpack('DRAG',desvar.ndim)
    force_x_model = Surfpack('FORCE_X',desvar.ndim)
    force_z_model = Surfpack('FORCE_Z',desvar.ndim)
    moment_y_model = Surfpack('MOMENT_Y',desvar.ndim)

    if (max_dsn == "ENRICHED"):

        # Load Data
        lift_model.load_data(os.path.join(project_folder,'enriched_points_lift.dat'))
        drag_model.load_data(os.path.join(project_folder,'enriched_points_drag.dat'))
        force_x_model.load_data(os.path.join(project_folder,'enriched_points_force_x.dat'))
        force_z_model.load_data(os.path.join(project_folder,'enriched_points_force_z.dat'))
        moment_y_model.load_data(os.path.join(project_folder,'enriched_points_moment_y.dat'))

    elif (max_dsn == "BUILD"):

        # Load Data
        lift_model.load_data(os.path.join(project_folder,'build_points_lift.dat'))
        drag_model.load_data(os.path.join(project_folder,'build_points_drag.dat'))
        force_x_model.load_data(os.path.join(project_folder,'build_points_force_x.dat'))
        force_z_model.load_data(os.path.join(project_folder,'build_points_force_z.dat'))
        moment_y_model.load_data(os.path.join(project_folder,'build_points_moment_y.dat'))

    else:

        project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
        project.compile_designs()

        if (max_dsn == "ALL"):
            dsn_number = len(project.designs)
        else:
            dsn_number = int(max_dsn)

        count = 0
        for index in range(dsn_number):
            design_container = project.designs[index]
            design = design_container.design
            if not design is None:
                dvs = desvar.pack(design.config)
                funcs = design.funcs
                if hasattr(funcs,'LIFT') and hasattr(funcs,'DRAG') and hasattr(funcs,'MOMENT_Y'):
                    count += 1
                    lift_model.add(dvs,funcs.LIFT)
                    drag_model.add(dvs,funcs.DRAG)
                    force_x_model.add(dvs,funcs.FORCE_X)
                    force_z_model.add(dvs,funcs.FORCE_Z)
                    moment_y_model.add(dvs,funcs.MOMENT_Y)
                    # print 'done ' + design.folder
                else:
                    print 'missing ' + design.folder
            else:
                print 'missing design'
        print count
        # Save Data
        lift_model.save_data(os.path.join(project_folder,'build_points_lift.dat'))
        drag_model.save_data(os.path.join(project_folder,'build_points_drag.dat'))
        force_x_model.save_data(os.path.join(project_folder,'build_points_force_x.dat'))
        force_z_model.save_data(os.path.join(project_folder,'build_points_force_z.dat'))
        moment_y_model.save_data(os.path.join(project_folder,'build_points_moment_y.dat'))

    # Build Model
    lift_model.build('kriging')
    drag_model.build('kriging')
    force_x_model.build('kriging')
    force_z_model.build('kriging')
    moment_y_model.build('kriging')

    # Save Model
    lift_model.save_model(os.path.join(project_folder,'model_lift.sps'))
    drag_model.save_model(os.path.join(project_folder,'model_drag.sps'))
    force_x_model.save_model(os.path.join(project_folder,'model_force_x.sps'))
    force_z_model.save_model(os.path.join(project_folder,'model_force_z.sps'))
    moment_y_model.save_model(os.path.join(project_folder,'model_moment_y.sps'))

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()

