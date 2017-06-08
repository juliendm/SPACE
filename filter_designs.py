#!/usr/bin/env python2.7

import os, time, sys, shutil, copy, glob
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
    parser.add_option("-a", "--apply", dest="apply", default="False",
                      help="apply", metavar="APPLY")
                      
    (options, args)=parser.parse_args()
    options.apply = options.apply.upper() == 'TRUE'

    filter_designs( options.project_folder, options.apply )

#: main()

def filter_designs( project_folder, apply = False ):

    ls = glob.glob(os.path.join(project_folder,'DESIGNS','*'))
    ls.sort()

    for dsn in ls:
        path = os.path.join(dsn,'fluid_surface_flow.dat')
        if not os.path.exists(path):
            print 'missing ' + dsn
            if apply: shutil.rmtree(dsn)

    ls_filtered = glob.glob(os.path.join(project_folder,'DESIGNS','*'))
    ls_filtered.sort()

    for index, dsn in enumerate(ls_filtered):
        new_index = '%03d' % (index+1)
        new_dsn = os.path.join(project_folder,'DESIGNS','DSN_'+new_index)
        if (not dsn == new_dsn) and apply:
            print 'moving ' + dsn
            shutil.move(dsn,new_dsn)

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()

