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

    project_folder = 'RESPONSE_SURFACE_DV_SUP'
    project = SPACE.io.load_data(project_folder + '/project.pkl')
    project.compile_designs()

    for design_container in project.designs:
        design = design_container.design

        if not design is None and  hasattr(design.funcs,'LIFT'):
            funcs = design.funcs
            print design.folder, funcs.LIFT, funcs.DRAG, funcs.MOMENT_Y
        else:
            print 'missing'

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()

