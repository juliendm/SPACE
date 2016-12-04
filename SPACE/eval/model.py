#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import time, os, gc, sys, shutil, copy, math
import numpy as np

from .. import surfpack   as spacesurfpack

# ----------------------------------------------------------------------
#  Model Class
# ----------------------------------------------------------------------

class Model(object):

    def __init__(self, config, state=None, folder='MODELS/MDL_*'):

        ndim = 5

    def eval(self):

        ndim = 5

    def enrich(self):

        ndim = 5        

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


