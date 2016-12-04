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
    project = SPACE.project.Project(config, state, folder='GEOMETRY_DVS_BIS')
    # project = SPACE.io.load_data('GEOMETRY_DVS_BIS/project.pkl')
    # print project.designs[0].design

    konfig = copy.deepcopy(config)

    procs = []

    # konfig.DV1 = '-0.5'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '-0.4'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '-0.3'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '-0.2'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '-0.1'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.1'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.2'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.3'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.4'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.5'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))

    # konfig.DV1 = '0.0'; konfig.DV2 = '-0.5'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '-0.4'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '-0.3'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '-0.2'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '-0.1'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.1'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.2'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.3'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.4'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.5'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))

    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.5'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.4'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.3'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.2'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.1'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.1'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.2'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.3'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.4'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.5'
    procs.append(project.func('STRUCTURE', konfig))

    for proc in procs:
        if not proc is None: proc.wait()

    project.deep_compile()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
