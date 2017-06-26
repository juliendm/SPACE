#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import time, os, gc, sys, shutil, copy, math
import numpy as np
import scipy as sp

import pyOpt

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

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_eval_sub:
            return self.lift_model_sub.eval(dvs)
        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_eval_sup:
            return self.lift_model_sup.eval(dvs)
        else:
            return self.interpolate_eval(dvs, self.lift_model_sub, self.lift_model_sup)

    def drag(self, dvs):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_eval_sub:
            return self.drag_model_sub.eval(dvs)
        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_eval_sup:
            return self.drag_model_sup.eval(dvs)
        else:
            return self.interpolate_eval(dvs, self.drag_model_sub, self.drag_model_sup)

    def force_z(self, dvs):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_eval_sub:
            return self.force_z_model_sub.eval(dvs)
        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_eval_sup:
            return self.force_z_model_sup.eval(dvs)
        else:
            return self.interpolate_eval(dvs, self.force_z_model_sub, self.force_z_model_sup)

    def force_z_raw_gradient(self, dvs):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_eval_sub:
            return self.force_z_model_sub.gradient(dvs)
        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_eval_sup:
            return self.force_z_model_sup.gradient(dvs)
        else:
            return self.interpolate_gradient(dvs,self.force_z_model_sub,self.force_z_model_sup)

    def moment_y(self, dvs, ref_origin_x):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_eval_sub:
            return self.moment_y_model_sub.eval(dvs) + (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * self.force_z_model_sub.eval(dvs)
        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_eval_sup:
            return self.moment_y_model_sup.eval(dvs) + (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * self.force_z_model_sup.eval(dvs)
        else:
            return self.interpolate_eval(dvs,self.moment_y_model_sub,self.moment_y_model_sup,ref_origin_x,self.force_z_model_sub,self.force_z_model_sup)

    def moment_y_raw_gradient(self, dvs, ref_origin_x):

        if dvs[self.desvar.mach_index] <= self.desvar.max_mach_eval_sub:
            return self.moment_y_model_sub.gradient(dvs) + (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * self.force_z_model_sub.gradient(dvs)
        elif dvs[self.desvar.mach_index] >= self.desvar.min_mach_eval_sup:
            return self.moment_y_model_sup.gradient(dvs) + (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * self.force_z_model_sup.gradient(dvs)
        else:
            return self.interpolate_gradient(dvs,self.moment_y_model_sub,self.moment_y_model_sup,ref_origin_x,self.force_z_model_sub,self.force_z_model_sup)








    def static_margin(self, dvs, ref_origin_x):

        moment_y_raw_grad = self.moment_y_raw_gradient(dvs,ref_origin_x)
        force_z_raw_grad = self.force_z_raw_gradient(dvs) 

        return - moment_y_raw_grad[self.desvar.aoa_index] / force_z_raw_grad[self.desvar.aoa_index] # * self.desvar.chain_rule_aoa(dvs)/self.desvar.chain_rule_aoa(dvs) CHANGE OF VARIABLES chain rule cancels out

    def k_alpha(self, dvs, ref_origin_x, static_margin_aug = 0.02):

        moment_y_raw_grad = self.moment_y_raw_gradient(dvs,ref_origin_x)
        force_z_raw_grad = self.force_z_raw_gradient(dvs)

        # eps = 0.00001
        # dvs_fd = copy.copy(dvs)

        # dvs_fd = self.desvar.forward_variable_change(dvs_fd)
        # dvs_fd[self.desvar.aoa_index] = dvs_fd[self.desvar.aoa_index] + eps
        # dvs_fd = self.desvar.reverse_variable_change(dvs_fd)
        # print '------------------'
        # print (self.force_z(dvs_fd)-self.force_z(dvs))/eps
        # print force_z_raw_grad[self.desvar.aoa_index]*self.desvar.chain_rule_aoa(dvs)

        # dvs_fd = self.desvar.forward_variable_change(dvs_fd)
        # dvs_fd[self.desvar.el_index] = dvs_fd[self.desvar.el_index] + eps
        # dvs_fd = self.desvar.reverse_variable_change(dvs_fd)
        # print '------------------'
        # print (self.moment_y(dvs_fd,ref_origin_x)-self.moment_y(dvs,ref_origin_x))/eps
        # print moment_y_raw_grad[self.desvar.el_index]*np.pi/180.0

        return (self.static_margin(dvs,ref_origin_x) - static_margin_aug) / ( (moment_y_raw_grad[self.desvar.el_index]*self.desvar.chain_rule_el) / (force_z_raw_grad[self.desvar.aoa_index]*self.desvar.chain_rule_aoa(dvs)) )

    def trim(self, dvs, ref_origin_x):

        new_dvs = copy.copy(dvs)

        def trim_function(x):
            self.deflections_update(new_dvs,x[0])
            return self.moment_y(new_dvs,ref_origin_x)

        x0 = sp.optimize.fsolve(trim_function, 0.0)
        self.deflections_update(new_dvs,x0[0])

        return new_dvs

    def deflections_update(self, dvs, x):

        if x < self.desvar.bf_bound[0]:
            dvs[self.desvar.bf_index] = self.desvar.bf_bound[0]
            dvs[self.desvar.el_index] = x - self.desvar.bf_bound[0]
        elif x > self.desvar.bf_bound[1]:
            dvs[self.desvar.bf_index] = self.desvar.bf_bound[1]
            dvs[self.desvar.el_index] = x - self.desvar.bf_bound[1]
        else:
            dvs[self.desvar.bf_index] = x
            dvs[self.desvar.el_index] = 0.0

    def max_trimmed_efficiency(self, dvs, ref_origin_x):

        new_dvs = copy.copy(dvs)

        def objfunc(x):

            new_dvs[self.desvar.bf_index] = x[0]
            new_dvs[self.desvar.el_index] = x[1]

            f = -self.lift(new_dvs)/self.drag(new_dvs) # minus sign to Maximize

            g = [self.moment_y(new_dvs,ref_origin_x)]

            fail = 0

            return f,g,fail

        opt_prob = pyOpt.Optimization('Max Trimmed Efficiency', objfunc)

        opt_prob.addObj('Efficiency')
        opt_prob.addCon('Trim', type='e', equal=0.0)

        opt_prob.addVar('bf','c',lower=self.desvar.bf_bound[0],upper=self.desvar.bf_bound[1],value=0.0)
        opt_prob.addVar('el','c',lower=self.desvar.el_bound[0],upper=self.desvar.el_bound[1],value=0.0)

        opt = pyOpt.SLSQP()
        opt.setOption('IPRINT',-1)
        opt.setOption('ACC',1e-5)

        [Y_min,X_min,Info] = opt(opt_prob, sens_type='FD') # Could Improve sens_type !!!!!!!!!!!

        return -Y_min[0]






    def interpolate_eval(self, dvs, model_sub, model_sup, ref_origin_x=0.0, model_sub_bis=None, model_sup_bis=None):

        dvs_sub = copy.copy(dvs)
        dvs_sub[self.desvar.mach_index] = self.desvar.max_mach_eval_sub
        dvs_sup = copy.copy(dvs)
        dvs_sup[self.desvar.mach_index] = self.desvar.min_mach_eval_sup

        coeff_sub = model_sub.eval(dvs_sub)
        coeff_sup = model_sup.eval(dvs_sup)

        if (model_sub_bis): coeff_sub += (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * model_sub_bis.eval(dvs_sub)
        if (model_sup_bis): coeff_sup += (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * model_sup_bis.eval(dvs_sup)

        return coeff_sub + (dvs[self.desvar.mach_index]-self.desvar.max_mach_eval_sub) * (coeff_sup-coeff_sub)/(self.desvar.min_mach_eval_sup-self.desvar.max_mach_eval_sub)

    def interpolate_gradient(self, dvs, model_sub, model_sup, ref_origin_x=0.0, model_sub_bis=None, model_sup_bis=None):

        dvs_sub = copy.copy(dvs)
        dvs_sub[self.desvar.mach_index] = self.desvar.max_mach_eval_sub
        dvs_sup = copy.copy(dvs)
        dvs_sup[self.desvar.mach_index] = self.desvar.min_mach_eval_sup

        coeff_sub = model_sub.gradient(dvs_sub)
        coeff_sup = model_sup.gradient(dvs_sup)

        if (model_sub_bis): coeff_sub += (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * model_sub_bis.gradient(dvs_sub)
        if (model_sup_bis): coeff_sup += (ref_origin_x-self.ref_origin_x_ini)/self.ref_length_moment * model_sup_bis.gradient(dvs_sup)

        return coeff_sub + (dvs[self.desvar.mach_index]-self.desvar.max_mach_eval_sub) * (coeff_sup-coeff_sub)/(self.desvar.min_mach_eval_sup-self.desvar.max_mach_eval_sub)


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


