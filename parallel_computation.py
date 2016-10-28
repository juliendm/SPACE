#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Config
    config = SPACE.io.Config('config.cfg')
    # State
    state  = SPACE.io.State()
    #state.find_files(config)

    # Project
    project = SPACE.project.Project(config, state, None, 'GEOMETRY_DVS')

    konfig = copy.deepcopy(config)
    project.func('STRUCTURE', konfig)

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
