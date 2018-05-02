#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np
import pwd
import pickle

from optparse import OptionParser

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack
from SPACE.util import LHC_unif, DesignVariables

sys.path.append(os.environ['SU2_RUN'])
import SU2

def main():

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--project", dest="project_folder",
                      help="project folder", metavar="PROJECT_FOLDER")
    parser.add_option("-r", "--regime", dest="regime", default="ON",
                      help="regime", metavar="REGIME")
    parser.add_option("-i", "--initiate", dest="initiate", default="False",
                      help="initiate", metavar="INITIATE")
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")


    (options, args)=parser.parse_args()
    options.initiate = options.initiate.upper() == 'TRUE'
    options.partitions = int( options.partitions )

    response_surface( options.filename       ,
                      options.project_folder ,
                      options.regime         ,
                      options.initiate       ,
                      options.partitions  )

#: main()

def response_surface( filename          ,
                      project_folder    ,
                      regime = 'BOTH'   ,
                      initiate = False  ,
                      partitions  = 0  ):

    # Project

    if os.path.exists(project_folder):
        project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
        project.compile_designs()
        config = project.config
    else:
        config = SPACE.io.Config(filename)
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    print '%d design(s) so far' % len(project.designs)

    config.NUMBER_PART = 0

    # Design Variables
    desvar = DesignVariables()

    # Load All Aero Models
    #load_aero_models()


    if regime == 'ON':
        XB = desvar.XB_STRUCT_ON
        ndim_struct = desvar.ndim_struct_on
        unpack_structure = desvar.unpack_structure_on
        pack_structure = desvar.pack_structure_on

        dry_mass_index = 3
        fuel_mass_index = 4
        thrust_index = 5
        pdyn_index = 6

        dv1_index = 7
        dv2_index = 8
        dv3_index = 9
        dv4_index = 10
        dv5_index = 11
        dv6_index = 12

        dvs_baseline = [1.1,0.0,1.0,20.0e3,20.0e3,2.5e6,20.0e3,0.0,0.0,0.0,0.0,0.0,0.0]

    elif regime == 'OFF':

        XB = desvar.XB_STRUCT_OFF
        ndim_struct = desvar.ndim_struct_off
        unpack_structure = desvar.unpack_structure_off
        pack_structure = desvar.pack_structure_off

    elif regime == 'BOTH':

        XB = desvar.XB_STRUCT
        ndim_struct = desvar.ndim_struct
        unpack_structure = desvar.unpack_structure
        pack_structure = desvar.pack_structure

        mach_index = 0
        rey_index = 1
        aoa_index = 2

        nx_index = 3
        nz_index = 4

        thrust_index = 5
        pdyn_index = 6
        fuel_mass_index = 7

        dv1_index = 8
        dv2_index = 9
        dv3_index = 10
        dv4_index = 11
        dv5_index = 12
        dv6_index = 13

        dvs_baseline = [1.1,0.0,1.0,3.0,-3.0,1.5e6,20.0e3,20.0e3,0.0,0.0,0.0,0.0,0.0,0.0]

        # Midpoint
        dvs_baseline = [4.5,0.0,0.0,3.0,-3.0,1.5e6,20.0e3,14.0e3,0.0,0.0,0.0,0.0,0.0,0.0]


    if initiate:


        vals = np.linspace(-0.5, 0.5, 10.0)
        new_dvs_vec = []
        for val in vals:
            dvs = copy.copy(dvs_baseline)
            dvs[dv2_index] = val
            new_dvs_vec.append(dvs)


        procs = []

        for new_dvs in new_dvs_vec:

            print 'New dvs: ', new_dvs

            konfig = copy.deepcopy(config)
            unpack_structure(konfig, new_dvs)

            proc = project.func('STRUCTURE', konfig)
            procs.append(proc)

        for proc in procs:
            proc.wait()


        vals = np.linspace(-0.5, 0.5, 10.0)
        new_dvs_vec = []
        for val in vals:
            dvs = copy.copy(dvs_baseline)
            dvs[dv3_index] = val
            new_dvs_vec.append(dvs)


        procs = []

        for new_dvs in new_dvs_vec:

            print 'New dvs: ', new_dvs

            konfig = copy.deepcopy(config)
            unpack_structure(konfig, new_dvs)

            proc = project.func('STRUCTURE', konfig)
            procs.append(proc)

        for proc in procs:
            proc.wait()


        # # nd = ndim_struct * 10
        # nd = 300

        # dvs_struct_filename = 'dvs_struct_' + regime + '.dat'

        # # X = LHC_unif(XB,nd)
        # # np.savetxt(dvs_struct_filename,X)

        # X = np.loadtxt(dvs_struct_filename)

        # for index_1 in range(nd/partitions):
        #     procs = []
        #     for index_2 in range(partitions):
        #         index = index_1*partitions+index_2

        #         dvs = X[index]

        #         konfig = copy.deepcopy(config)
        #         unpack_structure(konfig, dvs)

        #         proc = project.func('STRUCTURE', konfig)

        #         # force_redo_dsn_folder = None #'DSN_001'
        #         # proc = project.func('STRUCTURE', konfig, force_redo_dsn_folder)

        #         procs.append(proc)


        #     for proc in procs:
        #         proc.wait()

    else:

        na = 20
        number_per_ite = 12

#        flag = 'DRY_MASS'     # WHICH ONE ???????
        flag = 'STRUCTURE_MASS'

        threshold = 1.0001



        build_points_folder = os.path.join(project_folder,'BUILD_POINTS')

        model = Surfpack(flag, ndim_struct)
        model.load_data(os.path.join(build_points_folder,'build_points_' + flag + '.dat'))
        exclude = []

        next_dvs_vec_filename = 'next_dvs_vec_' + regime + '.dat'




        # new_dvs_vec = []
        # for dsn_index in range(12):
        #     design_folder = os.path.join(project_folder,'DESIGNS/DSN_%03d' % (dsn_index+1))
        #     local_config = SPACE.io.Config(os.path.join(design_folder,'config_DSN.cfg'))
        #     local_dvs = pack_structure(local_config)
        #     new_dvs_vec.append(local_dvs)
        # np.savetxt(next_dvs_vec_filename, new_dvs_vec)


        print 'Begin Iterating:'

        for ite in range(na):

            print 'Ite:', ite

            # NEW POINTS

            if ite == 0:

                new_dvs_vec = np.loadtxt(next_dvs_vec_filename)

            else:

                model.build('kriging')
                print 'Model built'

                new_dvs_vec = model.max_variance(XB, number = number_per_ite, exclude = exclude)
                np.savetxt(next_dvs_vec_filename, new_dvs_vec)

            effective_number_per_ite = len(new_dvs_vec)

            # COMPUTE

            procs = []

            number_design_before = len(project.designs)

            for new_dvs in new_dvs_vec:

                print 'New dvs: ', new_dvs

                konfig = copy.deepcopy(config)
                unpack_structure(konfig, new_dvs)

                proc = project.func('STRUCTURE', konfig)
                procs.append(proc)

            for proc in procs:
                proc.wait()

            number_design_after = len(project.designs)

            # READ RESULTS


#           for dsn_index in range(number_design_after-effective_number_per_ite,number_design_after):

            model = Surfpack(flag + '_' + str(ite), ndim_struct)
#            model.load_data(os.path.join(build_points_folder,'build_points_' + flag + '.dat'))
            exclude = []
            for dsn_index in range(0,number_design_after):

                design_folder = os.path.join(project_folder,'DESIGNS/DSN_%03d' % (dsn_index+1))

                history_file_name = os.path.join(design_folder,'STRUCTURE/history_structure.dat')
                postpro_file_name = os.path.join(design_folder,'STRUCTURE/postpro_load_1.dat')
                mass_file_name = os.path.join(design_folder,'STRUCTURE/lc0_mass_member.dat') # TO CHECK IF DONE INDEED

                if os.path.exists(history_file_name) and os.path.exists(postpro_file_name) and os.path.exists(mass_file_name):

                    history = np.loadtxt(history_file_name,delimiter=',',skiprows=1)

                    half_structure_mass = history[-1,1]

                    local_config = SPACE.io.Config(os.path.join(design_folder,'config_DSN.cfg'))
                    local_dvs = pack_structure(local_config)

                    # ADD VALUE

                    check = history[-1,2:5]
                    if (check[0] < threshold) and (check[1] < threshold) and (check[2] < threshold):

                        if flag == 'STRUCTURE_MASS':
                            model.add(local_dvs, half_structure_mass)
                        elif flag == 'DRY_MASS':
                            raise ValueError('Do not use, not correctly defined')

                    else:
                        exclude.append(local_dvs)
                        print 'Warning:', dsn_index+1
                else:
                    print 'Missing:', dsn_index+1


            model.save_data(os.path.join(build_points_folder,'enriched_points_' + flag + '.dat'))


        # Save Model

        model.build('kriging')
        model.save_model(os.path.join(project_folder,'model_' + flag + '.sps'))


#: response_surface()




def load_aero_models():

    n_models = 2180

    models_folder = 'RESPONSE_SURFACE_DV_SUP/DESIGNS/MODELS'

    for index in range(n_models):
        cp_model = Surfpack('CP_%05d' % (index+1), desvar.ndim)
        cp_model.load_model(os.path.join(models_folder,'model_cp_%05d.sps' % (index+1)))

    for index in range(n_models):
        cfx_model = Surfpack('CFX_%05d' % (index+1), desvar.ndim)
        cfx_model.load_model(os.path.join(models_folder,'model_cfx_%05d.sps' % (index+1)))

    for index in range(n_models):
        cfy_model = Surfpack('CFY_%05d' % (index+1), desvar.ndim)
        cfy_model.load_model(os.path.join(models_folder,'model_cfy_%05d.sps' % (index+1)))

    for index in range(n_models):
        cfz_model = Surfpack('CFZ_%05d' % (index+1), desvar.ndim)
        cfz_model.load_model(os.path.join(models_folder,'model_cfz_%05d.sps' % (index+1)))


#: load_models()







            #     vals = np.linspace(-0.5, 0.5, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[dv1_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 1:

            #     vals = np.linspace(0.0, 26.0e3, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[fuel_mass_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 2:

            #     vals = np.linspace(0.1e6, 3.0e6, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[thrust_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 3:

            #     vals = np.linspace(1.0e-6, 35.0e3, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[pdyn_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 4:

            #     vals = np.linspace(1.1, 8.0, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[mach_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 5:

            #     vals = np.linspace(-1.0, 1.0, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[aoa_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 6:

            #     vals = np.linspace(-3.0, 6.0, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[nx_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 7:

            #     vals = np.linspace(-6.0, 3.0, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[nz_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 8:

            #     vals = np.linspace(-0.5, 0.5, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[dv2_index] = val
            #         new_dvs_vec.append(dvs)

            # elif ite == 9:

            #     vals = np.linspace(-0.5, 0.5, 10.0)
            #     new_dvs_vec = []
            #     for val in vals:
            #         dvs = copy.copy(dvs_baseline)
            #         dvs[dv3_index] = val
            #         new_dvs_vec.append(dvs)








# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
