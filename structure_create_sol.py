#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, glob, re
# from .. import io   as spaceio
# from .  import func as spacefunc
# from ..io import redirect_folder, save_data

import subprocess
import numpy

SPACE_RUN = os.environ['SPACE_RUN']

from optparse import OptionParser

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE import io   as spaceio
from SPACE.eval import func as spacefunc
from SPACE.eval import design as spacedesign
from SPACE.io import redirect_folder, save_data

from SPACE.util import LHC_unif, DesignVariables
from SPACE.surfpack import Surfpack

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from sklearn.naive_bayes import GaussianNB

# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------

def create_sol(regime = 'ON'):

    n_models = 132
    models_folder = 'MODELS'
    build_points_folder = 'BUILD_POINTS'
    
    desvar = DesignVariables()


    if regime == 'ON':

        ndim_struct = desvar.ndim_struct_on
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

    elif regime == 'OFF':

        ndim_struct = desvar.ndim_struct_off
        pack_structure = desvar.pack_structure_off

        dry_mass_index = 3
        pdyn_index = 4

        dv1_index = 5
        dv2_index = 8
        dv3_index = 7
        dv4_index = 8
        dv5_index = 9
        dv6_index = 10

    elif regime == 'BOTH':

        ndim_struct = desvar.ndim_struct
        pack_structure = desvar.pack_structure

        thrust_index = 5
        pdyn_index = 6

        dv1_index = 8
        dv2_index = 9
        dv3_index = 10
        dv4_index = 11
        dv5_index = 12
        dv6_index = 13


    # mass_models = []
    # for index in range(n_models):
    #     mass_model = Surfpack('MASS_%05d' % (index+1), ndim_struct)
    #     mass_model.load_model(os.path.join(models_folder,'model_mass_%05d.sps' % (index+1)))
    #     mass_models.append(mass_model)

    # area_models = []
    # for index in range(n_models):
    #     area_model = Surfpack('AREA_%05d' % (index+1), ndim_struct)
    #     area_model.load_model(os.path.join(models_folder,'model_area_%05d.sps' % (index+1)))
    #     area_models.append(area_model)

    # thickness_models = []
    # thickness_clfs = []
    # thickness_values = []
    # for index in range(n_models):

    #     dat_file = os.path.join(build_points_folder,'enriched_points_thickness_%05d.dat' % (index+1))
    #     clf, values = get_clf(dat_file)
    #     thickness_clfs.append(clf)
    #     thickness_values.append(values)

    #     filename = os.path.join(models_folder,'model_thickness_filtered_%05d.sps' % (index+1))
    #     if os.path.exists(filename):
    #         thickness_model = Surfpack('THICKNESS_%05d' % (index+1), ndim_struct)
    #         thickness_model.load_model(filename)
    #         thickness_models.append(thickness_model)
    #     else:
    #         thickness_models.append(None)

    structure_mass_model = Surfpack('STRUCTURE_MASS', ndim_struct)
    structure_mass_model.load_model(os.path.join(models_folder,'model_structure_mass.sps'))

    structure_area_model = Surfpack('STRUCTURE_AREA', ndim_struct)
    structure_area_model.load_model(os.path.join(models_folder,'model_structure_area.sps'))

    structure_mass_model_538 = Surfpack('STRUCTURE_MASS_538', ndim_struct)
    structure_mass_model_538.load_model(os.path.join(models_folder,'model_structure_mass_538.sps'))


    # if regime == 'ON':
    #     print '///////////////////////////////////'
    #     dvs = [8,  -8.154462e-01 , -3.769461e-01,   2.105852e+04 ,  2.353362e+04  , 1.514796e+06 ,  2.065276e+04  , 1.156752e-01 , -1.911024e-01,  -4.267155e-01,  -1.732892e-01 , -3.821894e-01 ,  6.623392e-02]
    #     dvs = [6.681148e+00 ,  3.875544e-02 , -4.776890e-01 ,  5.181418e+03  , 1.119053e+04 ,  1.802179e+05 ,  6.625015e+03 ,  4.811014e-01 , -1.691460e-01 , -1.736698e-01 ,  2.441538e-01 , -3.857268e-01 ,  7.727754e-03]
    #     print structure_mass_model.eval(dvs)
    #     half_structure_mass_compounded = 0
    #     for index in range(n_models):
    #         area = area_models[index].eval(dvs)
    #         thickness = thickness_models[index].eval(dvs)
    #         if thickness < 0.0016: thickness = 0.0016
    #         half_structure_mass_compounded += area*thickness*2780.0
    #     print half_structure_mass_compounded
    #     dvs = [3 ,  1.000000e+00,  -4.428224e-01 ,  3.000000e+04  , 0.000000e+00  , 1.000000e+05 ,  1.000000e-06 , -5.000000e-01,  -5.000000e-01 , -5.000000e-01 , -2.000000e-01 , -5.000000e-01 , -2.000000e-01]
    #     print structure_mass_model.eval(dvs)
    #     half_structure_mass_compounded = 0
    #     for index in range(n_models):
    #         half_structure_mass_compounded += mass_models[index].eval(dvs)
    #     print half_structure_mass_compounded
    #     dvs = [8 ,  4.778873e-01   ,7.637604e-01 ,  6.573100e+03 ,  1.946402e+03  , 2.520920e+06 ,  2.220040e+04 , -3.427057e-01 ,  9.163581e-02 , -1.100301e-01  , 3.888815e-01 ,  4.278542e-01 ,  3.149009e-01]
    #     print structure_mass_model.eval(dvs)
    #     half_structure_mass_compounded = 0
    #     for index in range(n_models):
    #         half_structure_mass_compounded += mass_models[index].eval(dvs)
    #     print half_structure_mass_compounded 
    #     print '///////////////////////////////////'


#    for val in numpy.linspace(-0.5,0.5,10.0):
#        dvs[dv2_index] = val
#        mass = 0
#        for index in range(n_models):
#            mass += mass_models[index].eval(dvs)
#        print mass
#        print mass_models[49].eval(dvs)



    if regime == 'ON':
        #dvs = [1.1,0.0,1.0,20.0e3,20.0e3,2.5e6,20.0e3,0.0,0.0,0.0,0.0,0.0,0.0]
        dvs_baseline = [1.1,0.0,1.0,20.0e3,20.0e3,2.5e6,20.0e3,0.0,0.0,0.0,0.0,0.0,0.0]
        #dvs_baseline = [6.681148e+00 ,  3.875544e-02 , -4.776890e-01 ,  5.181418e+03  , 1.119053e+04 ,  1.802179e+05 ,  6.625015e+03 ,  4.811014e-01 , -1.691460e-01 , -1.736698e-01 ,  2.441538e-01 , -3.857268e-01 ,  7.727754e-03]
    elif regime == 'OFF':
        #dvs = [3.5,0.0,0.0,20.0e3,20.0e3,0.0,0.0,0.0,0.0,0.0,0.0]
        dvs_baseline = [1.1,0.0,1.0,20.0e3,20.0e3,0.0,0.0,0.0,0.0,0.0,0.0]
    elif regime == 'BOTH':
        #dvs = [3.5,0.0,0.0,20.0e3,20.0e3,0.0,0.0,0.0,0.0,0.0,0.0]
        dvs_baseline = [1.1,0.0,1.0,3.0,-3.0,1.5e6,20.0e3,20.0e3,0.0,0.0,0.0,0.0,0.0,0.0]


    # print 'FIG 1'

    # fig = plt.figure()
    # vals = numpy.linspace(0.1e6,3.0e6, 10.0)
    # for dv_index in [thrust_index]:
    #     half_structure_masses = []
    #     half_structure_masses_compounded = []
    #     for val in vals:
    #         dvs = copy.copy(dvs_baseline)
    #         dvs[dv_index] = val
    #         #dvs[dry_mass_index] = 20.e3
    #         half_structure_mass = structure_mass_model.eval(dvs)
    #         #effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass, dvs, regime)
    #         # for index in range(10):
    #         #     dvs[dry_mass_index] = 2.0*effective_half_dry_mass
    #         #     half_structure_mass = structure_mass_model.eval(dvs)
    #         #     effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass, dvs, regime)
    #         half_structure_masses.append(half_structure_mass)

    #         #dvs[dry_mass_index] = 20.e3

    #         # half_structure_mass_compounded = 0
    #         # for index in range(n_models):
    #         #     area = area_models[index].eval(dvs)
    #         #     thickness = thickness_models[index].eval(dvs)
    #         #     if thickness < 0.0016: thickness = 0.0016
    #         #     half_structure_mass_compounded += area*thickness*2780.0

    #         if thickness_models[51]: half_structure_mass_compounded = thickness_models[51].eval(dvs)

    #         #effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass_compounded, dvs, regime)
    #         # for index in range(10):
    #         #     dvs[dry_mass_index] = 2.0*effective_half_dry_mass
    #         #     half_structure_mass_compounded = 0
    #         #     for index in range(n_models):
    #         #         half_structure_mass_compounded += mass_models[index].eval(dvs)
    #         #     effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass_compounded, dvs, regime)
    #         half_structure_masses_compounded.append(half_structure_mass_compounded)

    #     #plt.plot(vals, half_structure_masses) #, color='red')
    #     plt.plot(vals, half_structure_masses_compounded) #, color='red')

    # plt.legend(['thrust'])
    # fig.savefig('fig_1.png')
    # plt.close(fig)

    # print 'FIG 2'

    # fig = plt.figure()
    # vals = numpy.linspace(0.0, 35.0e3, 10.0)
    # for dv_index in [pdyn_index]:
    #     half_structure_masses = []
    #     half_structure_masses_compounded = []
    #     for val in vals:

    #         dvs = copy.copy(dvs_baseline)

    #         dvs[dv_index] = val

    #         #dvs[dry_mass_index] = 20.e3
    #         half_structure_mass = structure_mass_model.eval(dvs)
    #         #effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass, dvs, regime)
    #         # for index in range(10):
    #         #     dvs[dry_mass_index] = 2.0*effective_half_dry_mass
    #         #     half_structure_mass = structure_mass_model.eval(dvs)
    #         #     effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass, dvs, regime)
    #         half_structure_masses.append(half_structure_mass)

    #         #dvs[dry_mass_index] = 20.e3
    #         half_structure_mass_compounded = 0
    #         for index in range(n_models):
    #             area = area_models[index].eval(dvs)
    #             thickness = thickness_models[index].eval(dvs)
    #             if thickness < 0.0016: thickness = 0.0016
    #             half_structure_mass_compounded += area*thickness*2780.0
    #         #effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass_compounded, dvs, regime)
    #         # for index in range(10):
    #         #     dvs[dry_mass_index] = 2.0*effective_half_dry_mass
    #         #     half_structure_mass_compounded = 0
    #         #     for index in range(n_models):
    #         #         half_structure_mass_compounded += mass_models[index].eval(dvs)
    #         #     effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass_compounded, dvs, regime)
    #         half_structure_masses_compounded.append(half_structure_mass_compounded)


    #     plt.plot(vals, half_structure_masses) #, color='red')
    #     plt.plot(vals, half_structure_masses_compounded) #, color='red')
    # plt.legend(['pdyn'])
    # fig.savefig('fig_2.png')
    # plt.close(fig)

    print 'FIG 3'

    fig = plt.figure()
    vals = numpy.linspace(-0.5, 0.5, 10.0)



    # eval_indexes = []
    # dvs = copy.copy(dvs_baseline)
    # for index in range(n_models):
    #     proba = thickness_clfs[index].predict_proba([dvs])
    #     thickness = thickness_values[index][numpy.argmin(proba)]
    #     if thickness == 0.0016: eval_indexes.append(index)  


    # eval_indexes = range(132) #eval_indexes[0:3]
    
    # print len(eval_indexes)
    # print eval_indexes

    # dat_file = os.path.join(build_points_folder,'enriched_points_mass_00001.dat')
    # data = numpy.loadtxt(dat_file)
    # intermediate_masses = [0.0 for index in range(len(data))]
    # for index in eval_indexes:
    #     dat_file = os.path.join(build_points_folder,'enriched_points_mass_%05d.dat' % (index+1))
    #     data = numpy.loadtxt(dat_file)
    #     for index_data in range(len(intermediate_masses)):
    #         val = data[index_data][-1]
    #         intermediate_masses[index_data] += val

    # intermediate_mass_model = Surfpack('INTERMEDIATE_MASS', ndim_struct)
    # dat_file = os.path.join(build_points_folder,'enriched_points_mass_00001.dat')
    # data = numpy.loadtxt(dat_file)
    # for index_data in range(0,len(data)):
    #     dvs = data[index_data][0:-1]
    #     intermediate_mass_model.add(dvs, intermediate_masses[index_data])

    # print 'BUILDING'
    # intermediate_mass_model.build('kriging')


    # dat_file = os.path.join(build_points_folder,'enriched_points_mass_00001.dat')
    # data = numpy.loadtxt(dat_file)
    # dvs = data[0][0:-1]

    # eval_indexes = []
    # for index in range(n_models):
    #     proba = thickness_clfs[index].predict_proba([dvs])
    #     thickness = thickness_values[index][numpy.argmin(proba)]
    #     if thickness == 0.0016: eval_indexes.append(index)  

    # print '/////////////////////////////////////////////'
    # print dvs

    # intermediate_masses = 0.0
    # for index in eval_indexes:
    #     dat_file = os.path.join(build_points_folder,'enriched_points_area_%05d.dat' % (index+1))
    #     data = numpy.loadtxt(dat_file)
    #     val = data[0][-1]
    #     intermediate_masses += val
    # print intermediate_masses


    # half_structure_mass_compounded = 0.0
    # for index in eval_indexes:
    #     area = area_models[index].eval(dvs)
    #     proba = thickness_clfs[index].predict_proba([dvs])
    #     thickness = thickness_values[index][numpy.argmin(proba)]
    #     half_structure_mass_compounded += area #*thickness*2780.0
    # print half_structure_mass_compounded
    # print '/////////////////////////////////////////////'




    for dv_index in [dv1_index,dv2_index,dv3_index]:
        print dv_index
        half_structure_masses = []
        half_intermediate_masses = []
        half_structure_masses_compounded = []
        for val in vals:

            dvs = copy.copy(dvs_baseline)

            dvs[dv_index] = val

            half_structure_mass = structure_mass_model.eval(dvs)
            half_structure_masses.append(half_structure_mass)

            # half_intermediate_mass = intermediate_mass_model.eval(dvs)
            # half_intermediate_masses.append(half_intermediate_mass)

            # half_structure_mass_compounded = 0
            # for index in eval_indexes: # range(n_models):
            #     area = area_models[index].eval(dvs)

            #     proba = thickness_clfs[index].predict_proba([dvs])
            #     thickness = thickness_values[index][numpy.argmin(proba)]
            #     if thickness:
            #         half_structure_mass_compounded += area*thickness*2780.0
            # half_structure_masses_compounded.append(half_structure_mass_compounded)


        plt.plot(vals, half_structure_masses) #, color='red')
        #plt.plot(vals, half_intermediate_masses) #, color='red')
        #plt.plot(vals, half_structure_masses_compounded) #, color='red')
    plt.legend(['dv1','dv1','dv1','dv2','dv2','dv2','dv3','dv3','dv3'])


            #effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass_compounded, dvs, regime)
            # for index in range(10):
            #     dvs[dry_mass_index] = 2.0*effective_half_dry_mass
            #     half_structure_mass_compounded = 0
            #     for index in range(n_models):
            #         half_structure_mass_compounded += mass_models[index].eval(dvs)
            #     effective_half_dry_mass = compute_effective_half_dry_mass(half_structure_mass_compounded, dvs, regime)

    print 'Save'

    fig.savefig('fig_3.png')
    plt.close(fig)


def get_clf(dat_file):

    data = numpy.loadtxt(dat_file)

    X_train = []
    y_train = []
    for index_data in range(len(data)):
        dvs = data[index_data][0:-1]
        val = data[index_data][-1]

        X_train.append(dvs.tolist())
        if val >= 0.1:
            y_train.append(2)
        elif val > 0.0016:
            y_train.append(1)
        else:
            y_train.append(0)

    X_train = numpy.array(X_train)
    y_train = numpy.array(y_train)

    min_y_train = min(y_train)
    max_y_train = max(y_train)

    if min_y_train == 0 and max_y_train == 2:
        values = [0.0016,None,0.1]
    elif min_y_train == 0:
        values = [0.0016,None]
    else:
        values = [None,0.1]

    clf = GaussianNB()

    clf.fit(X_train, y_train)

    return clf, values

def compute_effective_half_dry_mass(half_structure_mass, dvs, regime):


    if regime == 'ON':

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

    elif regime == 'OFF':

        dry_mass_index = 3
        pdyn_index = 4

        dv1_index = 5
        dv2_index = 8
        dv3_index = 7
        dv4_index = 8
        dv5_index = 9
        dv6_index = 10

    elif regime == 'BOTH':

        thrust_index = 5
        pdyn_index = 6

        dv1_index = 8
        dv2_index = 9
        dv3_index = 10
        dv4_index = 11
        dv5_index = 12
        dv6_index = 13



    half_dry_mass_kg_current_step = 0.5*dvs[dry_mass_index]

    max_fuel_mass = 26000.0
    half_mass_payload = 0.5 * 5000.0


    if regime == 'ON':
        fuel_percentage = dvs[fuel_mass_index] / max_fuel_mass
        half_max_thrust_newtons =  0.5 * dvs[thrust_index]
    elif regime == 'OFF':
        fuel_percentage = 0.0
        half_max_thrust_newtons =  0.0

    traj_max_pdyn_inf = dvs[pdyn_index]


    half_mass_fuel_kero = 0.5 * max_fuel_mass * 0.4
    half_mass_fuel_lox = 0.5 * max_fuel_mass * 0.6

    pounds_to_kg = 0.453592
    newtons_to_pounds = 0.224809
    kg_to_pounds = 2.20462
    meters_to_feet = 3.28084

    weight_pounds_current_step = 2.0*(half_dry_mass_kg_current_step+half_mass_fuel_kero+half_mass_fuel_lox)*kg_to_pounds
    max_thrust_pounds = 2.0*half_max_thrust_newtons*newtons_to_pounds
    body_length_feet = 18.0*meters_to_feet
    lox_density = 1141 # kg/m3
    kero_density = 810 # kg/m3
    N_engines = 1
    rocket_expansion_ratio = 77.5
    sts_tank_radius = 4.2 # m
    sts_tank_height = 46.9 # m
    sts_tank_area = 2.0*numpy.pi*sts_tank_radius*sts_tank_height + 2.0*numpy.pi*sts_tank_radius**2.0 # m2 # Assume: Right Cylinder
    sts_tank_empty_mass = 26500.0 # kg
    tank_mass_per_area = sts_tank_empty_mass/sts_tank_area # kg/m2
    tank_mass_per_area *= 0.6 # COMPOSITE REDUCE MASS BY 40 %
    technology_improvement = 0.2 # KEEP ONLY 20 % OF MASS

    # Landing Gear Weight

    weight_gear = 0.00916*weight_pounds_current_step**1.124
    half_mass_gear = weight_gear*0.5*pounds_to_kg

    # Avionics Weight

    weight_avionics = 66.37*weight_pounds_current_step**0.361
    half_mass_avionics = weight_avionics*0.5*pounds_to_kg * technology_improvement

    # Electrical System Weight

    weight_elec = 1.167*weight_pounds_current_step**0.5*body_length_feet**0.25
    half_mass_elec = weight_elec*0.5*pounds_to_kg * technology_improvement

    # Equipment Weight

    weight_equip = 1000.0 + 0.01*weight_pounds_current_step ############## Maybe reduce fixed value
    half_mass_equip = weight_equip*0.5*pounds_to_kg * technology_improvement

    # Tank LOX Weight

    volume_lox = 2.0*half_mass_fuel_lox/lox_density # m3
    lox_tank_radius = 1.65 # m
    lox_tank_height = volume_lox/numpy.pi/lox_tank_radius/lox_tank_radius
    lox_tank_area = 2.0*numpy.pi*lox_tank_radius*lox_tank_height + 2.0*numpy.pi*lox_tank_radius**2.0 # m2 # Assume: Right Cylinder
    half_mass_tank_lox = 0.5*tank_mass_per_area*lox_tank_area # kg

    # Tank KERO Weight

    volume_kero = 2.0*half_mass_fuel_kero/kero_density # m3
    kero_tank_radius = 1.65 # m
    kero_tank_height = volume_kero/numpy.pi/kero_tank_radius/kero_tank_radius
    kero_tank_area = 2.0*numpy.pi*kero_tank_radius*kero_tank_height + 2.0*numpy.pi*kero_tank_radius**2.0 # m2 # Assume: Right Cylinder
    half_mass_tank_kero = 0.5*tank_mass_per_area*kero_tank_area # kg

    # Engine Weight

    weight_engine = 0.00766*max_thrust_pounds + 0.00033*max_thrust_pounds*rocket_expansion_ratio**0.5 + 130.0*N_engines
    half_mass_engine = weight_engine*0.5*pounds_to_kg

    # Surface Dependant
    half_mass_tps = 100.0               # TODO: RESPONSE SURFACE !!!!!!!!!!!!!!!!!!!!!
    half_mass_hydraulic = 100.0




    half_additional_mass = 0.0

    half_additional_mass += half_mass_gear
    half_additional_mass += half_mass_avionics
    half_additional_mass += half_mass_elec
    half_additional_mass += half_mass_equip
    half_additional_mass += half_mass_tank_lox
    half_additional_mass += half_mass_tank_kero
    half_additional_mass += half_mass_engine
    half_additional_mass += fuel_percentage*half_mass_fuel_lox
    half_additional_mass += fuel_percentage*half_mass_fuel_kero
    half_additional_mass += half_mass_payload
    half_additional_mass += half_mass_hydraulic
    half_additional_mass += half_mass_tps





    

    half_dry_mass = half_structure_mass + half_additional_mass - fuel_percentage*(half_mass_fuel_kero+half_mass_fuel_lox)

    return half_dry_mass


    # half_wet_mass = half_structure_mass + half_additional_mass + (1.0-fuel_percentage)*(half_mass_fuel_kero+half_mass_fuel_lox)
    # half_wet_mass = half_dry_mass + half_mass_fuel_kero + half_mass_fuel_lox






#    sol = open('interp_sol.sol', 'w')
#    sol.write('MeshVersionFormatted 2\n\nDimension 3\n\nSolAtVertices\n' + str(n_models) + '\n4 1 1 1 1\n\n')
#    for index in range(n_models):
#        sol.write('%s %s %s %s\n' % (cp_models[index].eval(dvs), cfx_models[index].eval(dvs), cfy_models[index].eval(dvs), cfz_models[index].eval(dvs)))
#    sol.write('\nEnd\n')
#    sol.close()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-r", "--regime", dest="regime", default="ON",
                      help="regime", metavar="REGIME")


    (options, args)=parser.parse_args()

    create_sol( options.regime )



