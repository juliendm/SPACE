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

    # Project
    project_folder = 'GEOMETRY_DVS_MPI'
    if os.path.exists(project_folder):
        project = SPACE.io.load_data(project_folder + '/project.pkl')
        project.compile_designs()
        config = project.config
    else:
        config = SPACE.io.Config('config.cfg')
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    print '%d design(s) so far' % len(project.designs)
    #print project.designs[0].design.config.DV1

    konfig = copy.deepcopy(config)

    procs = []

    konfig.DV1 = '-0.5'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '-0.4'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '-0.3'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '-0.2'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '-0.1'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.1'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.2'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.3'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.4'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.5'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))

    konfig.DV1 = '0.0'; konfig.DV2 = '-0.5'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '-0.4'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '-0.3'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '-0.2'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '-0.1'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.1'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.2'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.3'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.4'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))
    konfig.DV1 = '0.0'; konfig.DV2 = '0.5'; konfig.DV3 = '0.0'
    procs.append(project.func('STRUCTURE', konfig))

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
        if not proc is None: proc.Disconnect()

    project.compile_designs()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
