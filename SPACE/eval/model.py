#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import time, os, gc, sys, shutil, copy, math
import numpy as np
import scipy as sp

from .. import surfpack   as spacesurfpack
from .. import util       as spacesutil

# ----------------------------------------------------------------------
#  AeroModel Class
# ----------------------------------------------------------------------

class AeroModel(object):

    def __init__(self, config, state=None, folder='MODELS/MDL_*'):

        self.config = config

        self.desvar = spacesutil.DesignVariables()

        self.lift_model_sub = spacesurfpack.Surfpack('LIFT_SUB',self.desvar.ndim)
        self.lift_model_sub.load_model(self.config.MODEL_LIFT_SUB)
        self.lift_model_sup = spacesurfpack.Surfpack('LIFT_SUP',self.desvar.ndim)
        self.lift_model_sup.load_model(self.config.MODEL_LIFT_SUP)

        self.drag_model_sub = spacesurfpack.Surfpack('DRAG_SUB',self.desvar.ndim)
        self.drag_model_sub.load_model(self.config.MODEL_DRAG_SUB)
        self.drag_model_sup = spacesurfpack.Surfpack('DRAG_SUP',self.desvar.ndim)
        self.drag_model_sup.load_model(self.config.MODEL_DRAG_SUP)

        self.force_z_model_sub = spacesurfpack.Surfpack('FORCE_Z_SUB',self.desvar.ndim)
        self.force_z_model_sub.load_model(self.config.MODEL_FORCE_Z_SUB)
        self.force_z_model_sup = spacesurfpack.Surfpack('FORCE_Z_SUP',self.desvar.ndim)
        self.force_z_model_sup.load_model(self.config.MODEL_FORCE_Z_SUP)

        self.moment_y_model_sub = spacesurfpack.Surfpack('MOMENT_Y_SUB',self.desvar.ndim)
        self.moment_y_model_sub.load_model(self.config.MODEL_MOMENT_Y_SUB)
        self.moment_y_model_sup = spacesurfpack.Surfpack('MOMENT_Y_SUP',self.desvar.ndim)
        self.moment_y_model_sup.load_model(self.config.MODEL_MOMENT_Y_SUP)

        self.ref_origin_x_ini = float(self.config.REF_ORIGIN_MOMENT_X)
        self.ref_length_moment = float(self.config.REF_LENGTH_MOMENT)

    def lift(self, dvs):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_sub:

            return self.lift_model_sub.eval(dvs)

        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_sup:

            return self.lift_model_sup.eval(dvs)

        else:

            dvs_sub = copy.copy(dvs)
            dvs_sub[self.desvar.mach_index] = self.desvar.max_mach_sub
            dvs_sup = copy.copy(dvs)
            dvs_sup[self.desvar.mach_index] = self.desvar.min_mach_sup

            cl_sub = self.lift_model_sub.eval(dvs_sub)
            cl_sup = self.lift_model_sup.eval(dvs_sup)

            return cl_sub + (dvs[self.desvar.mach_index]-self.desvar.max_mach_sub) * (cl_sup-cl_sub)/(self.desvar.min_mach_sup-self.desvar.max_mach_sub);

    def drag(self, dvs):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_sub:

            return self.drag_model_sub.eval(dvs)
            
        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_sup:

            return self.drag_model_sup.eval(dvs)

        else:

            dvs_sub = copy.copy(dvs)
            dvs_sub[self.desvar.mach_index] = self.desvar.max_mach_sub
            dvs_sup = copy.copy(dvs)
            dvs_sup[self.desvar.mach_index] = self.desvar.min_mach_sup

            cd_sub = self.drag_model_sub.eval(dvs_sub)
            cd_sup = self.drag_model_sup.eval(dvs_sup)

            return cd_sub + (dvs[self.desvar.mach_index]-self.desvar.max_mach_sub) * (cd_sup-cd_sub)/(self.desvar.min_mach_sup-self.desvar.max_mach_sub);

    def force_z(self, dvs):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_sub:

            return self.force_z_model_sub.eval(dvs)
            
        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_sup:

            return self.force_z_model_sup.eval(dvs)

        else:

            dvs_sub = copy.copy(dvs)
            dvs_sub[self.desvar.mach_index] = self.desvar.max_mach_sub
            dvs_sup = copy.copy(dvs)
            dvs_sup[self.desvar.mach_index] = self.desvar.min_mach_sup

            cz_sub = self.force_z_model_sub.eval(dvs_sub)
            cz_sup = self.force_z_model_sup.eval(dvs_sup)

            return cz_sub + (dvs[self.desvar.mach_index]-self.desvar.max_mach_sub) * (cz_sup-cz_sub)/(self.desvar.min_mach_sup-self.desvar.max_mach_sub);


    def moment_y(self, dvs, ref_origin_x):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_sub:

            return self.moment_y_model_sub.eval(dvs) + (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * self.force_z_model_sub.eval(dvs)

        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_sup:

            return self.moment_y_model_sup.eval(dvs) + (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * self.force_z_model_sup.eval(dvs)

        else:

            dvs_sub = copy.copy(dvs)
            dvs_sub[self.desvar.mach_index] = self.desvar.max_mach_sub
            dvs_sup = copy.copy(dvs)
            dvs_sup[self.desvar.mach_index] = self.desvar.min_mach_sup

            cmy_sub = self.moment_y_model_sub.eval(dvs_sub) + (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * self.force_z_model_sub.eval(dvs_sub)
            cmy_sup = self.moment_y_model_sup.eval(dvs_sup) + (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * self.force_z_model_sup.eval(dvs_sub)

            return cmy_sub + (dvs[self.desvar.mach_index]-self.desvar.max_mach_sub) * (cmy_sup-cmy_sub)/(self.desvar.min_mach_sup-self.desvar.max_mach_sub);

    def static_margin(self, dvs, ref_origin_x):

        eps = 0.001

        config_dvs = self.desvar.forward_variable_change(dvs)
        config_dvs[self.desvar.aoa_index] = config_dvs[self.desvar.aoa_index] + eps
        fd_dvs = self.desvar.reverse_variable_change(config_dvs)

        return - (self.moment_y(fd_dvs,ref_origin_x) - self.moment_y(dvs,ref_origin_x)) / (self.force_z(fd_dvs) - self.force_z(dvs))

    def trim(self, dvs, ref_origin_x):

        new_dvs = copy.copy(dvs)

        def trim_function(x):
            new_dvs[self.desvar.trim_index] = x[0]
            return self.moment_y(new_dvs,ref_origin_x)

        x0 = sp.optimize.fsolve(trim_function, 0.0)

        new_dvs[self.desvar.trim_index] = x0[0]

        return new_dvs

# ----------------------------------------------------------------------
#  RangedModel Class
# ----------------------------------------------------------------------

class RangedModel(object):

    def __init__(self, data):

        self.sub_models = []
        self.sub_ranges = []

        for key in data.keys():
            ranges = data[key].ranges
            self.sub_ranges.append(ranges)
            ndim = len(ranges)
            model = spacesurfpack.Surfpack(key, ndim)
            model.load_model(data[key].sps)
            self.sub_models.append(model)

    def dump_aero_database(self, filename, mach_vec, aoa_vec, dv1, dv2, dv3):

        dat_file = open(filename, 'w')

        for mach in mach_vec:
            for aoa in aoa_vec:        
                if (mach < 1.0 and aoa > 15.0) or (mach > 1.0 and aoa < 0.0):
                    dat_file.write('%.16E %.16E                    NaN\n' % (mach,aoa))
                else:
                    model = self.find_model([mach, aoa, dv1, dv2, dv3]) # VERY SLOW
                    dat_file.write('%.16E %.16E %.16E\n' % (mach,aoa,model.eval([mach, aoa, dv1, dv2, dv3])))

        dat_file.close()

    def find_model(self, dvs):

        correct_model = None

        for index_model in xrange(len(self.sub_models)):
            correct_model = index_model
            for index_dv in xrange(len(dvs)):
                if not (self.sub_ranges[index_model][index_dv][0] <= dvs[index_dv] <= self.sub_ranges[index_model][index_dv][1]):
                    correct_model = None
                    break
            if not correct_model is None:
                break

        if not correct_model is None:
            return self.sub_models[correct_model]
        else:
            return None




    # def dump_aero(self, dv1, dv2, dv3):

    #     aero = open('inputs/mk1_data/soar/aedb_dvs/soar_dvs_a.aero', 'w')

    #     aero.write('\nCL0   MACH   AOA\n\n')
    #     for mach in self.mach_vec: aero.write(str(mach)+' ')
    #     aero.write('\n')
    #     for aoa in self.aoa_vec: aero.write(str(aoa)+' ')
    #     aero.write('\n\n')
    #     for aoa in self.aoa_vec:
    #         for mach in self.mach_vec:
    #             if mach < 1.0:
    #                 model = self.lift_sub_model
    #             else:
    #                 model = self.lift_sup_model
    #             if (mach < 1.0 and aoa > 15.0) or (mach > 1.0 and aoa < 0.0):
    #                 aero.write('999999 ')
    #             else:
    #                 aero.write(str(model.eval([mach, aoa, dv1, dv2, dv3]))+' ')
    #         aero.write('\n')
    #     aero.write('\nCD0   MACH   AOA\n\n')
    #     for mach in self.mach_vec: aero.write(str(mach)+' ')
    #     aero.write('\n')
    #     for aoa in self.aoa_vec: aero.write(str(aoa)+' ')
    #     aero.write('\n\n')
    #     for aoa in self.aoa_vec:
    #         for mach in self.mach_vec:
    #             if mach < 1.0:
    #                 model = self.drag_sub_model
    #             else:
    #                 model = self.drag_sup_model
    #             if (mach < 1.0 and aoa > 15.0) or (mach > 1.0 and aoa < 0.0):
    #                 aero.write('999999 ')
    #             else:
    #                 aero.write(str(model.eval([mach, aoa, dv1, dv2, dv3]))+' ')
    #         aero.write('\n')
    #     aero.write('\n')

    #     aero.close()


