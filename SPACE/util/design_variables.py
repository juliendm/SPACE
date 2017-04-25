#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np

class DesignVariables(object):

    def __init__(self, regime='SUP'):

        if regime.upper() == 'SUB':
            mach_range = [0.3, 0.95]
        else:
            mach_range = [1.1, 8.0]

        self.XB = np.array([

            mach_range,         # dv_mach
            [-1.0, 1.0],       # dv_re
            [-1.0, 1.0],       # dv_aoa

            [-0.25, 0.5],      # dv_bf
            [-0.35, 0.5],      # dv_el

            [-0.5, 0.5],       # dv_geo1
            [-0.5, 0.5],       # dv_geo2
            [-0.5, 0.5],       # dv_geo3
            [-0.2, 0.5],       # dv_geo4
            [-0.5, 0.5],       # dv_geo5
            [-0.2, 0.5]        # dv_geo6

        ])

        self.ndim = len(self.XB)

    def unpack(self, config, dvs):

        dv_mach = dvs[0]
        dv_re = dvs[1]
        dv_aoa = dvs[2]
        dv_bf = dvs[3]
        dv_el = dvs[4]
        dv_geo1 = dvs[5]
        dv_geo2 = dvs[6]
        dv_geo3 = dvs[7]
        dv_geo4 = dvs[8]
        dv_geo5 = dvs[9]
        dv_geo6 = dvs[10]

        Reynolds, AoA = self.forward_variable_change(dv_mach, dv_re, dv_aoa)

        config.DV1 = str(dv_geo1); config.DV2 = str(dv_geo2); config.DV3 = str(dv_geo3); config.DV4 = str(dv_geo4); config.DV5 = str(dv_geo5); config.DV6 = str(dv_geo6)
        config.ELEVON_DEF = str(dv_el*180.0/np.pi); config.BODY_FLAP_DEF = str(dv_bf*180.0/np.pi)
        config.MACH_NUMBER = str(dv_mach); config.REYNOLDS_NUMBER = str(Reynolds); config.AoA= str(AoA)

    def pack(self, config):

        dv_mach = float(config.MACH_NUMBER)
        dv_re, dv_aoa = self.reverse_variable_change(dv_mach, float(config.REYNOLDS_NUMBER), float(config.AoA))
        dv_bf = float(config.BODY_FLAP_DEF)*np.pi/180.0
        dv_el = float(config.ELEVON_DEF)*np.pi/180.0
        dv_geo1 = float(config.DV1)
        dv_geo2 = float(config.DV2)
        dv_geo3 = float(config.DV3)
        dv_geo4 = float(config.DV4)
        dv_geo5 = float(config.DV5)
        dv_geo6 = float(config.DV6)

        return [dv_mach,dv_re,dv_aoa,dv_bf,dv_el,dv_geo1,dv_geo2,dv_geo3,dv_geo4,dv_geo5,dv_geo6]

    def forward_variable_change(self, dv_mach, dv_re, dv_aoa):

        Reynolds = 10.0**(-3.0/8.0*dv_mach+7.0+dv_re)

        if dv_mach >= 1.0:
            AoA = 3.92857*dv_mach+3.57143+dv_aoa*(1.78571*dv_mach+5.71429) # Supersonic
        else:
            AoA = (dv_aoa+1.0)*7.5 # Subsonic

        return Reynolds, AoA

    def reverse_variable_change(self, dv_mach, Reynolds, AoA):

        dv_re = np.log10(Reynolds) + 3.0/8.0*dv_mach - 7.0
        
        if dv_mach >= 1.0:
            dv_aoa = (AoA - (3.92857*dv_mach+3.57143)) / (1.78571*dv_mach+5.71429) # Supersonic
        else:
            dv_aoa = AoA/7.5-1.0 # Subsonic

        return dv_re, dv_aoa

