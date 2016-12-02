#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy

from .. import io  as spaceio

from interface import MIS       as SPACE_MIS

MISSION_ANALYSIS_RUN = os.environ['MISSION_ANALYSIS_RUN']

# ----------------------------------------------------------------------
#  Mission Simulation
# ----------------------------------------------------------------------

def mission ( config, database ): 
    
    # local copy
    konfig = copy.deepcopy(config)

    
    dv1 = float(konfig.DV1)
    dv2 = float(konfig.DV2)
    dv3 = float(konfig.DV3)

    cfg = open('soar.cfg', 'w')
    cfg.write('%.11f %.11f %.11f\n' % (dv1,dv2,dv3))
    cfg.close()

    # output

    os.mkdir('output')
    os.mkdir('output/soar')
    os.mkdir('output/us')

    # integrator

    os.mkdir('inputs')

    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/trajectory_setup.txt', 'inputs/trajectory_setup.txt')

    shutil.copytree(MISSION_ANALYSIS_RUN + '/integrator/inputs/environment_data', 'inputs/environment_data')

    os.mkdir('inputs/scenario_data')
#        shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/scenario_data/gpops_scenario.txt', dsn_folder + '/inputs/scenario_data/gpops_scenario.txt')
    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/scenario_data/integration_parameters.txt', 'inputs/scenario_data/integration_parameters.txt')

    os.makedirs('inputs/mk1_data/soar/aedb_dvs')
    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/soar/aero_data_dvs.txt', 'inputs/mk1_data/soar/aero_data_dvs.txt')

#        shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/soar/mass_budget_v01.txt', dsn_folder + '/inputs/mk1_data/soar/mass_budget_v01.txt')
    mass_dry = (1736.0 + dv1*(-202.22) + dv2*326.64 + dv3*134.31) * 2.0 + ( 1500.0 + 50.0 + 750.0 ) * 2.0 # engine, tanks, gnc
    mass_us = 6500.0
    mass_lox_kero = 24100.0
    mass_budget_file = open('inputs/mk1_data/soar/mass_budget_v01.txt','w')
    mass_budget_file.write('\n&soar_mbudget\n\nm0         = %f,\nm_stg      = %f,\n\n/\n\n' % (mass_dry+mass_lox_kero+mass_us, mass_us))
    mass_budget_file.close()

    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/soar/propulsion_v01.txt', 'inputs/mk1_data/soar/propulsion_v01.txt')
    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/soar/aedb_dvs/Rcurv_for_traj_regular.dat', 'inputs/mk1_data/soar/aedb_dvs/Rcurv_for_traj_regular.dat')

    os.makedirs('inputs/mk1_data/us/aedb_01')
    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/us/aero_data_it_01.txt', 'inputs/mk1_data/us/aero_data_it_01.txt')
    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/us/mass_budget_v01.txt', 'inputs/mk1_data/us/mass_budget_v01.txt')
    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/us/propulsion_v02.txt', 'inputs/mk1_data/us/propulsion_v02.txt')
    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/us/aedb_01/aedb_01.aero', 'inputs/mk1_data/us/aedb_01/aedb_01.aero')
    shutil.copy(MISSION_ANALYSIS_RUN + '/integrator/inputs/mk1_data/us/aedb_01/rnose.dat', 'inputs/mk1_data/us/aedb_01/rnose.dat')

    # optimizer

    os.makedirs('externals/auxfiles/aero')
    shutil.copy(MISSION_ANALYSIS_RUN + '/optimizer/externals/auxfiles/aero/rnose.dat', 'externals/auxfiles/aero/rnose.dat')

    os.mkdir('externals/auxfiles/mass')
    mass_dry_file = open('externals/auxfiles/mass/mass_dry.dat','w')  
    mass_dry_file.write('%f\n' % mass_dry)
    mass_dry_file.close()

    os.mkdir('externals/auxfiles/spice')
    shutil.copy(MISSION_ANALYSIS_RUN + '/optimizer/externals/auxfiles/spice/naif0010.tls', 'externals/auxfiles/spice')
    shutil.copy(MISSION_ANALYSIS_RUN + '/optimizer/externals/auxfiles/spice/pck00010.tpc', 'externals/auxfiles/spice')

    # compute
    database.dump_dat(dv1, dv2, dv3)

    shutil.copy('externals/auxfiles/aero/CD0_a.dat', 'externals/auxfiles/aero/CD0_r.dat')
    shutil.copy('externals/auxfiles/aero/CL0_a.dat', 'externals/auxfiles/aero/CL0_r.dat')


    # Run Solution
    SPACE_MIS(konfig)


    # info out
    info = spaceio.State()
    # info.FUNCTIONS.update( aerodynamics )
    # info.FILES.DIRECT = konfig['RESTART_FLOW_FILENAME']
    # info.HISTORY.DIRECT = history
    
    return info
