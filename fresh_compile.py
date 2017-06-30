#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
from optparse import OptionParser
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-p", "--project", dest="project_folder",
                      help="project folder", metavar="PROJECT_FOLDER")
                      
    (options, args)=parser.parse_args()

    fresh_compile( options.project_folder )

#: main()

def fresh_compile( project_folder ):

    if os.path.exists(os.path.join(project_folder,'project.pkl')):
        project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
    else:
        config = SPACE.io.Config('config.cfg')
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    
    project.fresh_compile_designs(project_folder=project_folder)

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()

