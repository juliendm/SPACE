#!/usr/bin/env python2.7

import os, time, sys, shutil
import numpy as np

class DesignVariables(object):

    def __init__(self):

        self.ndim = 11

        self.mach_index = 0    # mach
        self.rey_index = 1     # reynolds
        self.aoa_index = 2     # aoa
        self.bf_index = 3      # bf
        self.el_index = 4      # el

        self.bf_bound = [-0.25, 0.5]
        self.el_bound = [-0.35, 0.5]

        self.XB_SUB = np.array([

            [0.3, 0.95],       # dv_mach
            [-1.0, 1.0],       # dv_rey
            [-1.0, 1.0],       # dv_aoa

            self.bf_bound,      # dv_bf
            self.el_bound,      # dv_el

            [-0.5, 0.5],       # dv_geo1
            [-0.5, 0.5],       # dv_geo2
            [-0.5, 0.5],       # dv_geo3
            [-0.2, 0.5],       # dv_geo4
            [-0.5, 0.5],       # dv_geo5
            [-0.2, 0.5]        # dv_geo6

        ])

        self.XB_SUP = np.array([

            [1.1, 8.0],       # dv_mach
            [-1.0, 1.0],       # dv_rey
            [-1.0, 1.0],       # dv_aoa

            self.bf_bound,      # dv_bf
            self.el_bound,      # dv_el

            [-0.5, 0.5],       # dv_geo1
            [-0.5, 0.5],       # dv_geo2
            [-0.5, 0.5],       # dv_geo3
            [-0.2, 0.5],       # dv_geo4
            [-0.5, 0.5],       # dv_geo5
            [-0.2, 0.5]        # dv_geo6

        ])


        self.ndim_struct_on = 13
        self.ndim_struct_off = 11

        self.XB_STRUCT_ON = np.array([
            self.XB_SUP[0],
            self.XB_SUP[1],
            self.XB_SUP[2],
            [ 5.0e3,  30.0e3],    # dv_dry_mass   [kg]
            [ 0.0,    26.0e3],       # dv_fuel_mass  [kg]    # WARNING: THIS HAS TO BE CONSISTENT WITH MAX_FUEL_MASS IN CONFIG FILE
            [ 0.1e6,   3.0e6],     # dv_thrust      [N]
            [ 0.0,    35.0e3],       # dv_pdyn       [Pa]
            self.XB_SUP[5],
            self.XB_SUP[6],
            self.XB_SUP[7],
            self.XB_SUP[8],
            self.XB_SUP[9],
            self.XB_SUP[10],
        ])

        self.XB_STRUCT_OFF = np.array([
            self.XB_STRUCT_ON[0],
            self.XB_STRUCT_ON[1],
            self.XB_STRUCT_ON[2],
            self.XB_STRUCT_ON[3],       # dv_dry_mass   [kg]
            self.XB_STRUCT_ON[6],       # dv_pdyn       [Pa]
            self.XB_STRUCT_ON[7],
            self.XB_STRUCT_ON[8],
            self.XB_STRUCT_ON[9],
            self.XB_STRUCT_ON[10],
            self.XB_STRUCT_ON[11],
            self.XB_STRUCT_ON[12],
        ])

        self.max_mach_eval_sub = 0.99
        self.min_mach_eval_sup = 1.01

        # to multiply with raw_grad
        self.chain_rule_mach = 1.0
        self.chain_rule_el = np.pi/180.0
        self.chain_rule_bf = np.pi/180.0

    def unpack(self, config, dvs):

        config_dvs = self.forward_variable_change(dvs)

        config.MACH_NUMBER = str(config_dvs[0]); config.REYNOLDS_NUMBER = str(config_dvs[1]); config.AoA= str(config_dvs[2])
        config.BODY_FLAP_DEF = str(config_dvs[3]); config.ELEVON_DEF = str(config_dvs[4])
        config.DV1 = str(config_dvs[5]); config.DV2 = str(config_dvs[6]); config.DV3 = str(config_dvs[7]); config.DV4 = str(config_dvs[8]); config.DV5 = str(config_dvs[9]); config.DV6 = str(config_dvs[10])

    def pack(self, config):

        return self.reverse_variable_change([float(config.MACH_NUMBER), float(config.REYNOLDS_NUMBER), float(config.AoA), float(config.BODY_FLAP_DEF), float(config.ELEVON_DEF), float(config.DV1),
        float(config.DV2), float(config.DV3), float(config.DV4), float(config.DV5), float(config.DV6)])

    def unpack_structure_on(self, config, dvs):

        config_dvs = self.forward_variable_change([dvs[0],dvs[1],dvs[2],0.0,0.0,dvs[7],dvs[8],dvs[9],dvs[10],dvs[11],dvs[12]])

        config.MACH_NUMBER = str(config_dvs[0]); config.REYNOLDS_NUMBER = str(config_dvs[1]); config.AoA= str(config_dvs[2])
        config.BODY_FLAP_DEF = str(config_dvs[3]); config.ELEVON_DEF = str(config_dvs[4])
        config.DV1 = str(config_dvs[5]); config.DV2 = str(config_dvs[6]); config.DV3 = str(config_dvs[7]); config.DV4 = str(config_dvs[8]); config.DV5 = str(config_dvs[9]); config.DV6 = str(config_dvs[10])

        config.DRY_MASS       = str(dvs[3]) # dv_dry_mass   [kg]
        config.FUEL_MASS      = str(dvs[4]) # dv_fuel_mass  [kg]
        config.THRUST         = str(dvs[5]) # dv_thrust      [N]
        config.P_DYN_INF      = str(dvs[6]) # dv_pdyn       [Pa]

    def pack_structure_on(self, config):

        dvs = self.reverse_variable_change([float(config.MACH_NUMBER), float(config.REYNOLDS_NUMBER), float(config.AoA), float(config.BODY_FLAP_DEF), float(config.ELEVON_DEF), float(config.DV1),
        float(config.DV2), float(config.DV3), float(config.DV4), float(config.DV5), float(config.DV6)])

        return [
            dvs[0],
            dvs[1],
            dvs[2],
            float(config.DRY_MASS),             # dv_dry_mass   [kg]
            float(config.FUEL_MASS),            # dv_fuel_mass  [kg]
            float(config.THRUST),               # dv_thrust      [N]
            float(config.P_DYN_INF),            # dv_pdyn       [Pa]
            dvs[5],
            dvs[6],
            dvs[7],
            dvs[8],
            dvs[9],
            dvs[10]
        ]

    def unpack_structure_off(self, config, dvs):

        config_dvs = self.forward_variable_change([dvs[0],dvs[1],dvs[2],0.0,0.0,dvs[5],dvs[6],dvs[7],dvs[8],dvs[9],dvs[10]])

        config.MACH_NUMBER = str(config_dvs[0]); config.REYNOLDS_NUMBER = str(config_dvs[1]); config.AoA= str(config_dvs[2])
        config.BODY_FLAP_DEF = str(config_dvs[3]); config.ELEVON_DEF = str(config_dvs[4])
        config.DV1 = str(config_dvs[5]); config.DV2 = str(config_dvs[6]); config.DV3 = str(config_dvs[7]); config.DV4 = str(config_dvs[8]); config.DV5 = str(config_dvs[9]); config.DV6 = str(config_dvs[10])

        config.DRY_MASS       = str(dvs[3]) # dv_dry_mass   [kg]
        config.FUEL_MASS      = 0.0         # dv_fuel_mass  [kg]
        config.THRUST         = 0.0         # dv_thrust      [N]
        config.P_DYN_INF      = str(dvs[4]) # dv_pdyn       [Pa]


    def pack_structure_off(self, config):

        dvs = self.reverse_variable_change([float(config.MACH_NUMBER), float(config.REYNOLDS_NUMBER), float(config.AoA), float(config.BODY_FLAP_DEF), float(config.ELEVON_DEF), float(config.DV1),
        float(config.DV2), float(config.DV3), float(config.DV4), float(config.DV5), float(config.DV6)])

        return [
            dvs[0],
            dvs[1],
            dvs[2],
            float(config.DRY_MASS),              # dv_dry_mass   [kg]
            float(config.P_DYN_INF),             # dv_pdyn       [Pa]
            dvs[5],
            dvs[6],
            dvs[7],
            dvs[8],
            dvs[9],
            dvs[10]
        ]

    def forward_variable_change(self, dvs):

        dv_mach = dvs[0]
        dv_rey = dvs[1]
        dv_aoa = dvs[2]
        dv_bf = dvs[3]
        dv_el = dvs[4]

        dv_geo1 = dvs[5]
        dv_geo2 = dvs[6]
        dv_geo3 = dvs[7]
        dv_geo4 = dvs[8]
        dv_geo5 = dvs[9]
        dv_geo6 = dvs[10]

        Reynolds = 10.0**(-3.0/8.0*dv_mach+7.0+dv_rey)

        if dv_mach >= 1.0:
            AoA = 3.92857*dv_mach+3.57143+dv_aoa*(1.78571*dv_mach+5.71429) # Supersonic
        else:
            AoA = (dv_aoa+1.0)*7.5 # Subsonic

        return [dv_mach, Reynolds, AoA, dv_bf*180.0/np.pi, dv_el*180.0/np.pi, dv_geo1, dv_geo2, dv_geo3, dv_geo4, dv_geo5, dv_geo6]

    def reverse_variable_change(self, config_dvs):

        dv_mach = config_dvs[0]
        Reynolds = config_dvs[1]
        AoA = config_dvs[2]
        Body_Flap = config_dvs[3]
        Elevon = config_dvs[4]

        dv_geo1 = config_dvs[5]
        dv_geo2 = config_dvs[6]
        dv_geo3 = config_dvs[7]
        dv_geo4 = config_dvs[8]
        dv_geo5 = config_dvs[9]
        dv_geo6 = config_dvs[10]

        dv_rey = np.log10(Reynolds) + 3.0/8.0*dv_mach - 7.0
        
        if dv_mach >= 1.0:
            dv_aoa = (AoA - (3.92857*dv_mach+3.57143)) / (1.78571*dv_mach+5.71429) # Supersonic
        else:
            dv_aoa = AoA/7.5-1.0 # Subsonic

        return [dv_mach, dv_rey, dv_aoa, Body_Flap*np.pi/180.0, Elevon*np.pi/180.0, dv_geo1, dv_geo2, dv_geo3, dv_geo4, dv_geo5, dv_geo6]

    def chain_rule_aoa(self, dvs):

        dv_mach = dvs[0]

        if dv_mach >= 1.0:
            return 1.0/(1.78571*dv_mach+5.71429) # Supersonic
        else:
            return 1.0/7.5 # Subsonic



