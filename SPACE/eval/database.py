#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import time, os, gc, sys, shutil, copy, math
import numpy as np

from .. import surfpack   as spacesurfpack

SURROGATE_COMPUTE = os.environ['SURROGATE_COMPUTE']

# ----------------------------------------------------------------------
#  Database Class
# ----------------------------------------------------------------------

class Database(object):

    def __init__(self, config, state=None, folder='DATABASE'):

        ndim = 5

        self.lift_sup_model = spacesurfpack.Surfpack('lift_sup', ndim)
        self.lift_sup_model.load_model(SURROGATE_COMPUTE + '/AERO/sps/sp_gp_model.lift_sup.sps')

        self.drag_sup_model = spacesurfpack.Surfpack('drag_sup', ndim)
        self.drag_sup_model.load_model(SURROGATE_COMPUTE + '/AERO/sps/sp_gp_model.drag_sup.sps')

        self.lift_sub_model = spacesurfpack.Surfpack('lift_sub', ndim)
        self.lift_sub_model.load_model(SURROGATE_COMPUTE + '/AERO/sps/sp_gp_model.lift_sub.sps')

        self.drag_sub_model = spacesurfpack.Surfpack('drag_sub', ndim)
        self.drag_sub_model.load_model(SURROGATE_COMPUTE + '/AERO/sps/sp_gp_model.drag_sub.sps')

        self.mach_vec = [0.5, 0.6, 0.7, 0.8, 0.9, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.5, 2.7, 3.0, 3.2, 3.5, 3.7, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0] # no more than 60 !!!!!!!!
        self.aoa_vec = [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0] # no more than 60 !!!!!!!!

    def dump_dat(self, dv1, dv2, dv3):

        lift = open('externals/auxfiles/aero/CL0_a.dat', 'w')

        for mach in self.mach_vec:
            if mach < 1.0:
                model = self.lift_sub_model
            else:
                model = self.lift_sup_model
            for aoa in self.aoa_vec:        
                if (mach < 1.0 and aoa > 15.0) or (mach > 1.0 and aoa < 0.0):
                    lift.write('%.16E %.16E                    NaN\n' % (mach,aoa))
                else:
                    lift.write('%.16E %.16E %.16E\n' % (mach,aoa,model.eval([mach, aoa, dv1, dv2, dv3])))


        lift.close()


        drag = open('externals/auxfiles/aero/CD0_a.dat', 'w')

        for mach in self.mach_vec:
            if mach < 1.0:
                model = self.drag_sub_model
            else:
                model = self.drag_sup_model
            for aoa in self.aoa_vec:
                if (mach < 1.0 and aoa > 15.0) or (mach > 1.0 and aoa < 0.0):
                    drag.write('%.16E %.16E                    NaN\n' % (mach,aoa))
                else:
                    drag.write('%.16E %.16E %.16E\n' % (mach,aoa,model.eval([mach, aoa, dv1, dv2, dv3])))
        drag.close()

    def dump_aero(self, dv1, dv2, dv3):

        aero = open('inputs/mk1_data/soar/aedb_dvs/soar_dvs_a.aero', 'w')

        aero.write('\nCL0   MACH   AOA\n\n')
        for mach in self.mach_vec: aero.write(str(mach)+' ')
        aero.write('\n')
        for aoa in self.aoa_vec: aero.write(str(aoa)+' ')
        aero.write('\n\n')
        for aoa in self.aoa_vec:
            for mach in self.mach_vec:
                if mach < 1.0:
                    model = self.lift_sub_model
                else:
                    model = self.lift_sup_model
                if (mach < 1.0 and aoa > 15.0) or (mach > 1.0 and aoa < 0.0):
                    aero.write('999999 ')
                else:
                    aero.write(str(model.eval([mach, aoa, dv1, dv2, dv3]))+' ')
            aero.write('\n')
        aero.write('\nCD0   MACH   AOA\n\n')
        for mach in self.mach_vec: aero.write(str(mach)+' ')
        aero.write('\n')
        for aoa in self.aoa_vec: aero.write(str(aoa)+' ')
        aero.write('\n\n')
        for aoa in self.aoa_vec:
            for mach in self.mach_vec:
                if mach < 1.0:
                    model = self.drag_sub_model
                else:
                    model = self.drag_sup_model
                if (mach < 1.0 and aoa > 15.0) or (mach > 1.0 and aoa < 0.0):
                    aero.write('999999 ')
                else:
                    aero.write(str(model.eval([mach, aoa, dv1, dv2, dv3]))+' ')
            aero.write('\n')
        aero.write('\n')

        aero.close()


