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
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
                      
    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )



    # Project
    project_folder = 'RESPONSE_SURFACE_DV'
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
    konfig.NUMBER_PART = options.partitions

    procs = []

    konfig.DV1 = '-0.5'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'; konfig.DV4 = '0.0'; konfig.DV5 = '0.0'; konfig.DV6 = '0.0'
    konfig.ELEVON_DEF = '0.0'; konfig.BODY_FLAP_DEF = '0.0'
    konfig.MACH_NUMBER= '0.8'; konfig.REYNOLDS_NUMBER = '10E6'; konfig.AoA= '10.0'
    proc = project.func('AERODYNAMICS', konfig)
    proc.wait()

    konfig.DV1 = '-0.25'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'; konfig.DV4 = '0.0'; konfig.DV5 = '0.0'; konfig.DV6 = '0.0'
    konfig.ELEVON_DEF = '0.0'; konfig.BODY_FLAP_DEF = '0.0'
    konfig.MACH_NUMBER= '3.0'; konfig.REYNOLDS_NUMBER = '8E5'; konfig.AoA= '20.0'
    proc = project.func('AERODYNAMICS', konfig)
    proc.wait()

    konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'; konfig.DV4 = '0.0'; konfig.DV5 = '0.0'; konfig.DV6 = '0.0'
    konfig.ELEVON_DEF = '0.0'; konfig.BODY_FLAP_DEF = '0.0'
    konfig.MACH_NUMBER= '3.0'; konfig.REYNOLDS_NUMBER = '6E5'; konfig.AoA= '20.0'
    proc = project.func('AERODYNAMICS', konfig)
    proc.wait()

    konfig.DV1 = '0.25'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'; konfig.DV4 = '0.0'; konfig.DV5 = '0.0'; konfig.DV6 = '0.0'
    konfig.ELEVON_DEF = '0.0'; konfig.BODY_FLAP_DEF = '0.0'
    konfig.MACH_NUMBER= '3.0'; konfig.REYNOLDS_NUMBER = '4E5'; konfig.AoA= '20.0'
    proc = project.func('AERODYNAMICS', konfig)
    proc.wait()

    konfig.DV1 = '0.5'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'; konfig.DV4 = '0.0'; konfig.DV5 = '0.0'; konfig.DV6 = '0.0'
    konfig.ELEVON_DEF = '0.0'; konfig.BODY_FLAP_DEF = '0.0'
    konfig.MACH_NUMBER= '3.0'; konfig.REYNOLDS_NUMBER = '2E5'; konfig.AoA= '20.0'
    proc = project.func('AERODYNAMICS', konfig)
    proc.wait()

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

    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.5'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.4'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.3'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.2'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '-0.1'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.0'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.1'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.2'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.3'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.4'
    # procs.append(project.func('STRUCTURE', konfig))
    # konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = '0.5'
    # procs.append(project.func('STRUCTURE', konfig))

    for proc in procs:
        #if not proc is None: proc.Disconnect()
        if not proc is None: proc.wait()

    project.compile_designs()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
