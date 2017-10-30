#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from __future__ import division

import os, sys, shutil, copy, time

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing, PGMbody, PGMshell
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone
from GeoMACH.PSM import Airframe

from scipy import interpolate
import numpy as np

from .. import io  as spaceio
from interface import MRG       as SPACE_MRG
from interface import SUR       as SPACE_SUR

SPACE_RUN = os.environ['SPACE_RUN']

# ----------------------------------------------------------------------
#  Geometry Simulation
# ----------------------------------------------------------------------

def geometry ( config ): 
    
    # local copy
    konfig = copy.deepcopy(config)

    # Run Solution
    pgm = Spaceplane()

    pgm.dv5 = float(konfig.DV5)
    pgm.dv6 = float(konfig.DV6)
    pgm.body_flap_deflection = float(konfig.BODY_FLAP_DEF)

    # FLUID

    pgm.structure = False
    pgm.tail = False

    wing_width_section_1_ini = 1.54
    wing_width_section_2_ini = 2.41
    pgm.wing_width_section_1, pgm.wing_width_section_2 = compute_wing_profiles(float(konfig.DV1),float(konfig.DV2),
        float(konfig.DV3),float(konfig.DV4),wing_width_section_1_ini,wing_width_section_2_ini,float(konfig.ELEVON_DEF),'fluid')

    bse = pgm.initialize()

    pgm.comps['lwing'].set_airfoil('profile_fluid_1.dat')
    pgm.comps['rwing'].set_airfoil('profile_fluid_37.dat')
    pgm.comps['flap'].set_airfoil('profile_cs.dat')
    pgm.compute_all()
    # pgm.compute_normals()
    # pgm.compute_all()

    # Export Mesh
    bse.vec['pt_str'].export_MESH(konfig.FLUID_SURFACE + '_back.mesh')
    SPACE_MRG(konfig)
    SPACE_SUR(konfig)

    # STRUCTURE
    
    if konfig.STRUCT != 'NONE':

        pgm.structure = True
        pgm.tail = True

        wing_width_section_1_ini = 1.54
        wing_width_section_2_ini = 2.41
        pgm.wing_width_section_1, pgm.wing_width_section_2 = compute_wing_profiles(float(konfig.DV1),float(konfig.DV2),
            float(konfig.DV3),float(konfig.DV4),wing_width_section_1_ini,wing_width_section_2_ini,float(konfig.ELEVON_DEF),'struct')

        bse = pgm.initialize()

        pgm.comps['lwing'].set_airfoil('profile_struct_1.dat')
        pgm.comps['rwing'].set_airfoil('profile_struct_37.dat')
        pgm.comps['flap'].set_airfoil('profile_cs.dat')
        pgm.compute_all()

        pgm.meshStructure(konfig.STRUCT)

    # info out
    info = spaceio.State()
    info.FILES.FLUID_SURFACE_MESH = konfig.FLUID_SURFACE + '.mesh'
    info.FILES.FLUID_SURFACE_BACK_MESH = konfig.FLUID_SURFACE + '_back.mesh'
    if konfig.STRUCT != 'NONE':
        info.FILES.STRUCT_BDF = konfig.STRUCT + '.bdf'
        info.FILES.STRUCT_MESH = konfig.STRUCT + '.mesh'
        info.FILES.STRUCT_SURFACE_MESH = konfig.STRUCT + '_surface.mesh'
    return info

class Spaceplane(PGMconfiguration):

    def _define_comps(this):

        if this.structure:
            this.fuse_num_x = 16
        else:
            this.fuse_num_x = 14

        this.comps['fuse'] = PGMbody(num_x=this.fuse_num_x, num_y=4, num_z=8) #num_x is the number of surfacees

        if this.structure:
            this.comps['lwing'] = PGMwing(num_x=11, num_z=12, left_closed=True)
            this.comps['rwing'] = PGMwing(num_x=11, num_z=12, right_closed=True)
        else:
            this.comps['lwing'] = PGMwing(num_x=9, num_z=12, left_closed=True)
            this.comps['rwing'] = PGMwing(num_x=9, num_z=12, right_closed=True)

        if this.tail: this.comps['ctail'] = PGMwing(num_x=3, num_z=4, left_closed=True)

        this.comps['flap'] = PGMwing(num_x=6, num_z=4, left_closed=True)

        this.comps['fuse_f'] = PGMcone(this, 'fuse', 'front', 7.0)
        this.comps['fuse_r'] = PGMcone(this, 'fuse', 'rear', 3.0)
        this.comps['lwing_t'] = PGMtip(this, 'lwing', 'left', 0.1)
        this.comps['rwing_t'] = PGMtip(this, 'rwing', 'right', 0.1)
        if this.tail: this.comps['ctail_t'] = PGMtip(this, 'ctail', 'left', 3.0)
      #  this.comps['flap_t'] = PGMtip(this, 'flap', 'left', 0.1)

        #this.comps['lwing_fuse'] = PGMjunction(this, 'fuse', 'bot', 'E', [0,1], 'lwing', 'right')
        if this.structure:
            this.comps['lwing_fuse'] = PGMjunction(this, 'fuse', 'lft', 'E', [3,2], 'lwing', 'right')
            this.comps['rwing_fuse'] = PGMjunction(this, 'fuse', 'rgt', 'W', [3,1], 'rwing', 'left')
        else:
            this.comps['lwing_fuse'] = PGMjunction(this, 'fuse', 'lft', 'E', [3,1], 'lwing', 'right')
            this.comps['rwing_fuse'] = PGMjunction(this, 'fuse', 'rgt', 'W', [3,2], 'rwing', 'left')

        if this.tail: this.comps['ctail_fuse'] = PGMjunction(this, 'fuse', 'top', 'E', [3,this.fuse_num_x-5], 'ctail', 'right')

        this.comps['flap_fuse'] = PGMjunction(this, 'fuse', 'bot', 'S', [this.fuse_num_x-2,0], 'flap', 'right')

    def _define_params(this):

        fuse = this.comps['fuse'].props #1st parameter: streamwise
        fuse['pos'].params[''] = PGMparameter(2, 3) #Specify the cylinder front and end
        fuse['nor'].params[''] = PGMparameter(1, 3) #Normal of each section
        fuse['scl'].params[''] = PGMparameter(1, 1)
        fuse['flt'].params[''] = PGMparameter(3, 4, pos_u=[0.0,0.2,1.0])
        fuse['pos'].params['nose'] = PGMparameter(3, 3, pos_u=[0,0.15,0.3], order_u=3)
        fuse['scl'].params['nose_y'] = PGMparameter(3, 3, pos_u=[0,0.15,0.3], order_u=3)
        fuse['scl'].params['nose_z'] = PGMparameter(3, 3, pos_u=[0,0.1,0.2], order_u=3)
        fuse['pos'].params['bottom'] = PGMparameter(3, 3, pos_u=[0.2,0.88,1.0], order_u=2)
        fuse['scl'].params['bottom'] = PGMparameter(3, 3, pos_u=[0.2,0.88,1.0], order_u=2)
        fuse['rot'].params['bottom'] = PGMparameter(3, 3, pos_u=[0.2,0.88,1.0], order_u=2)

        lwing = this.comps['lwing'].props
        lwing['pos'].params[''] = PGMparameter(1, 3)
        lwing['pos'].params['lin'] = PGMparameter(3, 3, order_u=2)
        lwing['scl'].params[''] = PGMparameter(3, 3, order_u=2)

        rwing = this.comps['rwing'].props
        rwing['pos'].params[''] = PGMparameter(1, 3)
        rwing['pos'].params['lin'] = PGMparameter(3, 3, order_u=2)
        rwing['scl'].params[''] = PGMparameter(3, 3, order_u=2)

        if this.tail:
            ctail = this.comps['ctail'].props
            ctail['pos'].params[''] = PGMparameter(1, 3)
            ctail['pos'].params['lin'] = PGMparameter(2, 3)
            ctail['scl'].params[''] = PGMparameter(2, 3)
            ctail['nor'].params[''] = PGMparameter(1, 3)

        flap = this.comps['flap'].props
        flap['pos'].params[''] = PGMparameter(1, 3)
        flap['pos'].params['lin'] = PGMparameter(2, 3)
        flap['scl'].params[''] = PGMparameter(2, 3)
        flap['rot'].params[''] = PGMparameter(2, 3)

    def _compute_params(this):

        nose_origin_x = -0.69

        # fuselage

        fuse_length = 18.0 # 18.0 #########################
        fuse_radius = 1.75
        fuse_bottom = 0.1

        fuse = this.comps['fuse'].props
        fuse['pos'].params[''].val([[0,0,0],[fuse_length+nose_origin_x,0,0]])
        fuse['nor'].params[''].val([1.0,0.0,0.0])
        fuse['scl'].params[''].val([fuse_radius])
        fuse['flt'].params[''].val([[0.3,0.3,0.1,0.1],[0.3,0.3,1.0,1.0],[0.3,0.3,1.0,1.0]])
        fuse['pos'].params['nose'].val([[0,-0.7,0],[0,0,0],[0,0,0]])
        fuse['scl'].params['nose_y'].val([[0,-1.0,0],[0,0,0],[0,0,0]])
        fuse['scl'].params['nose_z'].val([[-0.8,0,0],[0,0,0],[0,0,0]])
        #fuse['pos'].params['bottom'].val([[0,0,0],[0,-fuse_bottom,0],[0,-0.6,0]])
        fuse['pos'].params['bottom'].val([[0,0,0],[0,-fuse_bottom,0],[0,0.35,0]])
        #fuse['scl'].params['bottom'].val([[0,0,0],[0,fuse_bottom,0],[0,0.6,0]])
        fuse['scl'].params['bottom'].val([[0,0,0],[0,fuse_bottom,0],[0,-0.35,0]])
        fuse['rot'].params['bottom'].val([[0,0,0],[0,11,0],[0,11,0]])
        
        # wings

        wing_origin_x = 3.5
        wing_origin_y = -1.2
        wing_origin_z = 1.75

        wing_length = 20.0 # need to be 20 for consistency. Length of the wing is actually controlled in profiles.py

        lwing = this.comps['lwing'].props
        lwing['pos'].params[''].val([wing_origin_x+nose_origin_x,wing_origin_y,wing_origin_z])
        lwing['pos'].params['lin'].val([[0,0,0],[0,0,this.wing_width_section_1],[0,0,this.wing_width_section_1+this.wing_width_section_2]])
        lwing['scl'].params[''].val([[wing_length,1,1],[wing_length,1,1],[wing_length,1,1]])

        rwing = this.comps['rwing'].props
        rwing['pos'].params[''].val([wing_origin_x+nose_origin_x,wing_origin_y,-wing_origin_z-(this.wing_width_section_1+this.wing_width_section_2)])
        rwing['pos'].params['lin'].val([[0,0,0],[0,0,this.wing_width_section_2],[0,0,this.wing_width_section_1+this.wing_width_section_2]])
        rwing['scl'].params[''].val([[wing_length,1,1],[wing_length,1,1],[wing_length,1,1]])

        # tails

        if this.tail:
            tail_root_origin_x = 13.3 ################
            tail_tip_origin_x = 16.3 - 0.8 ################
            tail_length_root = 4.1 - 0.2 ################
            tail_length_tip = 2.5
            tail_height = 3.8

            ctail = this.comps['ctail'].props
            ctail['pos'].params[''].val([tail_root_origin_x+nose_origin_x,fuse_radius,0.0])
            ctail['pos'].params['lin'].val([[0,0,0],[tail_tip_origin_x-tail_root_origin_x,tail_height,0.0]])
            ctail['scl'].params[''].val([[tail_length_root,3.0,1.0],[tail_length_tip,3.0,1.0]])
            ctail['nor'].params[''].val([1.0,0.0,0.0])

        # flap

        flap_origin_x = 16.2 + this.dv6
        flap_origin_y = -1.8
        flap_width = 3.4
        flap_length = 2.0 + this.dv5

        flap_deflection = this.body_flap_deflection

        flap = this.comps['flap'].props
        flap['pos'].params[''].val([flap_origin_x+nose_origin_x,flap_origin_y,-flap_width/2.0])
        flap['pos'].params['lin'].val([[0,0,0],[flap_length * np.cos(-flap_deflection*np.pi/180.0), flap_length * np.sin(-flap_deflection*np.pi/180.0), 0]])
        flap['scl'].params[''].val([[flap_width,2.0,1],[flap_width,0.3,1]])
        flap['rot'].params[''].val([[180,90,0],[180,90,0]])

        return [], [], []

    def _set_bspline_options(this):
        comps = this.comps


        comps['fuse'].faces['bot'].set_option('num_pt', 'v', [8,10,15,10]+[6]*(this.fuse_num_x-4), both=False)

        comps['fuse'].faces['bot'].set_option('num_pt', 'u', [5,5,5,5,5,5,5,5], both=False)
        comps['fuse'].faces['lft'].set_option('num_pt', 'u', [6,6,6,8], both=False)

        comps['lwing'].faces['upp'].set_option('num_pt', 'v', [4,4,4,4,4,4,5,5,5,5,5,5], both=False)

        comps['flap'].faces['upp'].set_option('num_pt', 'v', [6,6,6,6], both=False)


    def meshStructure(this, filename):
        afm = Airframe(this, 1)

        
        middle_spar = 11

        hinge_spar = 5
        hinge_location = 0.22 # 0.285
        #idims = np.array([ 0.08, 0.15, 0.22, hinge_location, 0.3875,  0.49, 0.5925, 0.695, 0.7975, 0.9])
        # idims = np.linspace(0.08,0.9,9) #13
        idims = np.append(np.linspace(0.08,hinge_location,5)[0:-1],np.linspace(hinge_location,0.9,14))

        
        angle_rib = 3
        angle_location = 0.55 #20.0/37.0 #this.wing_width_section_1/this.wing_width_section_2 # + 0.05 #########
        jdims = np.array([0.0, 0.18333333, angle_location, 0.7, 0.85, 1.0])
        # jdims = np.append(np.linspace(0,angle_location,4)[0:-1],np.linspace(angle_location,1,4)) #13

        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                if i < hinge_spar and j < 1: continue
                afm.addVertFlip('MRIBW:%02d:l:%02d' % (j,i),'lwing',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('MRIBW:%02d:r:%02d' % (j,i),'rwing',[idims[i],1-jdims[j]],[idims[i+1],1-jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                if i < hinge_spar and j < 1: continue
                if i in [0, hinge_spar, idims.shape[0]-1]: # middle_spar, 
                    afm.addVertFlip('MSPARW:%02d:l:%02d' % (i,j),'lwing',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                    afm.addVertFlip('MSPARW:%02d:r:%02d' % (i,j),'rwing',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]])
                else:
                    if i in [6,8,10,12,14,16] and j >= angle_rib - 1: continue
                    afm.addVertFlip('MSTRINGW:%02d:a:l:%02d' % (i,j),'lwing',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('MSTRINGW:%02d:b:l:%02d' % (i,j),'lwing',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[0.15,0])
                    afm.addVertFlip('MSTRINGW:%02d:a:r:%02d' % (i,j),'rwing',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('MSTRINGW:%02d:b:r:%02d' % (i,j),'rwing',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]],w=[0.15,0])

        for i in range(hinge_spar,idims.shape[0]):
            if i in [hinge_spar, idims.shape[0]-1]:
                afm.addCtrVert('MSPARC:%02d' % (i),'lwing','rwing',idims[i])
            else:
                afm.addCtrVert('MSTRINGC:%02d:a' % (i),'lwing','rwing',idims[i],w=[1,0.85])
                afm.addCtrVert('MSTRINGC:%02d:b' % (i),'lwing','rwing',idims[i],w=[0.15,0])
        for i in range(hinge_spar,idims.shape[0]-1):
            afm.addCtr('MSKINC:a:%02d' % (i),'lwing','rwing',0,[idims[i],idims[i+1]])
        for i in range(hinge_spar,idims.shape[0]-1):
            afm.addCtr('MSKINC:b:%02d' % (i),'lwing','rwing',1,[1-idims[i],1-idims[i+1]])

        if this.tail:
            idims = np.array([0.2,0.4,0.6,0.8])
            jdims = np.linspace(0.0,0.95,10)
            for i in range(idims.shape[0]-1):
                for j in range(jdims.shape[0]):
                    afm.addVertFlip('MRIBV:%02d:%02d' % (j,i),'ctail',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
            for i in range(idims.shape[0]):
                for j in range(jdims.shape[0]-1):
                    afm.addVertFlip('MSPARV:%02d:%02d' % (i,j),'ctail',[idims[i],jdims[j]],[idims[i],jdims[j+1]])

        idims = np.linspace(0.15,0.85,4)
        jdims = np.linspace(0.0,1.0,8)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('MRIBF:%02d:%02d' % (j,i),'flap',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                afm.addVertFlip('MSPARF:%02d:%02d' % (i,j),'flap',[idims[i],jdims[j]],[idims[i],jdims[j+1]])

        idims = np.linspace(0,1,4) # 4
        jdims = np.linspace(0,1,20) # 13
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVert('MFRAME:%02d:1:%02d' % (j,i),'fuse',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.96],i=[0,2])
                afm.addVert('MFRAME:%02d:2:%02d' % (j,i),'fuse',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.96],i=[1,3])
                afm.addVert('MFRAME:%02d:3:%02d' % (j,i),'fuse',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.96],i=[2,0])
                afm.addVert('MFRAME:%02d:4:%02d' % (j,i),'fuse',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.96],i=[3,1])
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]-1):
                afm.addVert('MLONG:%02d:1:%02d' % (i,j),'fuse',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.96],i=[0,2])
                afm.addVert('MLONG:%02d:2:%02d' % (i,j),'fuse',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.96],i=[1,3])
                afm.addVert('MLONG:%02d:3:%02d' % (i,j),'fuse',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.96],i=[2,0])
                afm.addVert('MLONG:%02d:4:%02d' % (i,j),'fuse',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.96],i=[3,1])

        afm.preview('preview')
        afm.mesh()
        afm.computeMesh(filename)

def compute_wing_profiles(dv1, dv2, dv3, dv4, wing_width_section_1_ini, wing_width_section_2_ini, deflection, case = 'fluid'):

    n_profiles_0 = 4    
    n_profiles_1 = 16
    n_profiles_2 = 19

    scale_ini = 0.92 # 0.97 #

    x_hinge = 10.5 * scale_ini + dv4 # 9.78 * scale_ini # 
    y_hinge = -0.3

    profile = np.loadtxt(SPACE_RUN + '/SPACE/util/profiles/%s/edge_1.dat' % case)
    x_u_root = profile[0:65,0] * scale_ini
    x_l_root = profile[65:130,0] * scale_ini
    y_u_root = profile[0:65,1]
    y_l_root = profile[65:130,1]

    x_u_edge = np.zeros((3,65))
    y_u_edge = np.zeros((3,65))
    x_l_edge = np.zeros((3,65))
    y_l_edge = np.zeros((3,65))

    profile = np.loadtxt(SPACE_RUN + '/SPACE/util/profiles/%s/edge_3.dat' % case)
    x_u_edge[0,:] = profile[0:65,0] * scale_ini
    x_l_edge[0,:] = profile[65:130,0] * scale_ini
    y_u_edge[0,:] = profile[0:65,1]
    y_l_edge[0,:] = profile[65:130,1]

    profile = np.loadtxt(SPACE_RUN + '/SPACE/util/profiles/%s/edge_23.dat' % case)
    x_u_edge[1,:] = profile[0:65,0] * scale_ini
    x_l_edge[1,:] = profile[65:130,0] * scale_ini
    y_u_edge[1,:] = profile[0:65,1]
    y_l_edge[1,:] = profile[65:130,1]

    profile = np.loadtxt(SPACE_RUN + '/SPACE/util/profiles/%s/edge_45.dat' % case)
    x_u_edge[2,:] = profile[0:65,0] * scale_ini
    x_l_edge[2,:] = profile[65:130,0] * scale_ini
    y_u_edge[2,:] = profile[0:65,1]
    y_l_edge[2,:] = profile[65:130,1]

    profile = np.loadtxt(SPACE_RUN + '/SPACE/util/profiles/rae2822.dat')
    x_canvas = profile[0:65,0]





    # cutting_length = max(x_u_root)
    # cutting_index = 0
    # for k in range(len(x_u_edge[0,:])):
    #     if (x_u_edge[0,k] > cutting_length) :
    #         cutting_index = k
    #         break

    # a_u = (y_u_edge[0,cutting_index-1]-y_u_edge[0,cutting_index])/(x_u_edge[0,cutting_index-1]-x_u_edge[0,cutting_index])
    # b_u = y_u_edge[0,cutting_index] - a_u*x_u_edge[0,cutting_index]
    # x_u_in_cut_ini = x_u_edge[0,0:cutting_index+1]
    # x_u_in_cut_ini[-1] = cutting_length
    # y_u_in_cut_ini = y_u_edge[0,0:cutting_index+1]
    # y_u_in_cut_ini[-1] = a_u*x_u_in_cut_ini[-1] + b_u 

    # a_l = (y_l_edge[0,cutting_index-1]-y_l_edge[0,cutting_index])/(x_l_edge[0,cutting_index-1]-x_l_edge[0,cutting_index])
    # b_l = y_l_edge[0,cutting_index] - a_l*x_l_edge[0,cutting_index]
    # x_l_in_cut_ini = x_l_edge[0,0:cutting_index+1]
    # x_l_in_cut_ini[-1] = cutting_length
    # y_l_in_cut_ini = y_l_edge[0,0:cutting_index+1]
    # y_l_in_cut_ini[-1] = a_l*x_l_in_cut_ini[-1] + b_l

    # x_u_in_cut = x_canvas*(np.max(x_u_in_cut_ini)-np.min(x_u_in_cut_ini))+np.min(x_u_in_cut_ini)
    # x_l_in_cut = x_canvas*(np.max(x_l_in_cut_ini)-np.min(x_l_in_cut_ini))+np.min(x_l_in_cut_ini)
    # x_u_in_cut[0] = x_u_in_cut_ini[0]
    # x_l_in_cut[0] = x_l_in_cut_ini[0]
    # x_u_in_cut[-1] = x_u_in_cut_ini[-1]
    # x_l_in_cut[-1] = x_l_in_cut_ini[-1]

    # f_u = interpolate.interp1d(x_u_in_cut_ini,y_u_in_cut_ini)
    # f_l = interpolate.interp1d(x_l_in_cut_ini,y_l_in_cut_ini)
    # y_u_in_cut = f_u(x_u_in_cut)
    # y_l_in_cut = f_l(x_l_in_cut)

    # n_smooth = 12
    # y_target = y_l_in_cut[-1]
    # a_u_r = np.array([[x_u_in_cut[-n_smooth]*x_u_in_cut[-n_smooth],x_u_in_cut[-n_smooth],1],[x_u_in_cut[-1]*x_u_in_cut[-1],x_u_in_cut[-1],1],[2.0*x_u_in_cut[-n_smooth],1,0]])
    # b_u_r = np.array([y_u_in_cut[-n_smooth],y_target,(y_u_in_cut[-n_smooth]-y_u_in_cut[-n_smooth-1])/(x_u_in_cut[-n_smooth]-x_u_in_cut[-n_smooth-1])])
    # coeffs_u_r = np.linalg.solve(a_u_r,b_u_r)
    # for ind in range(n_smooth):
    #     index_r = -(ind+1)
    #     y_u_in_cut[index_r] = coeffs_u_r[0]*x_u_in_cut[index_r]*x_u_in_cut[index_r] + coeffs_u_r[1]*x_u_in_cut[index_r] + coeffs_u_r[2]


    x_u_in_cut = x_u_edge[0,:]
    y_u_in_cut = y_u_edge[0,:]

    x_l_in_cut = x_l_edge[0,:]
    y_l_in_cut = y_l_edge[0,:]



    # Design Variables

    a_strake = (min(x_u_edge[1,:])-min(x_u_edge[0,:]))/wing_width_section_1_ini
    b_strake = min(x_u_edge[1,:])-a_strake*wing_width_section_1_ini
    print a_strake, b_strake

    a_lead = (min(x_u_edge[2,:])-min(x_u_edge[1,:]))/wing_width_section_2_ini
    b_lead_ini = min(x_u_edge[1,:])-a_lead*wing_width_section_1_ini
    b_lead = b_lead_ini + dv1
    print a_lead, b_lead_ini

    a_trail = (max(x_u_edge[2,:])-max(x_u_edge[1,:]))/wing_width_section_2_ini
    b_trail_ini = max(x_u_edge[1,:])-a_trail*wing_width_section_1_ini
    b_trail = b_trail_ini + dv2 + dv4 # !!!!!!!!!!!!!!!!!!!!!!!!
    print a_trail, b_trail_ini

    wing_width_ini = wing_width_section_1_ini + wing_width_section_2_ini
    wing_width = wing_width_ini + dv3
    wing_width_section_1 = (b_lead-b_strake)/(a_strake-a_lead);
    wing_width_section_2 = wing_width - wing_width_section_1

    x_min = [0.0]*3
    x_min[0] = b_strake
    x_min[1] = a_lead*wing_width_section_1+b_lead
    x_min[2] = a_lead*wing_width+b_lead;

    x_max = [0.0]*3
    x_max[0] = b_trail
    x_max[1] = a_trail*wing_width_section_1+b_trail
    x_max[2] = a_trail*wing_width+b_trail

    for n in range(3):
        x_distrib = x_u_edge[n,:]-x_hinge
        pos_front_ini = min(x_distrib)
        pos_rear_ini = max(x_distrib)
        for k in range(len(x_distrib)):
            if (x_distrib[k] < 0):
                x_distrib[k] = x_distrib[k]/pos_front_ini*(x_min[n]-x_hinge)
            else:
                x_distrib[k] = x_distrib[k]/pos_rear_ini*(x_max[n]-x_hinge)

        x_u_edge[n,:] = x_distrib + x_hinge
        x_l_edge[n,:] = x_distrib + x_hinge


    # Deflection of Elevon

    angle = -deflection * np.pi / 180.0
    pc_hinge = 0.85

    filter_k = [0]*3
    filter_a = [0.0]*3
    filter_b = [0.0]*3

    for n in range(3):
        if angle >= 0.0:
            for i in range(len(x_u_edge[n,:])):
                if (x_u_edge[n,i] < x_hinge):
                    filter_k[n] = i
        else:
            for i in range(len(x_l_edge[n,:])):
                if (x_l_edge[n,i] < x_hinge):
                    filter_k[n] = i

    for n in range(3):
        if angle >= 0.0:
            index_filter = len(x_u_edge[n,:])-1
            filter_a[n] = (y_u_edge[n,index_filter] - y_u_edge[n,filter_k[n]]) / (x_u_edge[n,index_filter] - x_u_edge[n,filter_k[n]])
            filter_b[n] = y_u_edge[n,index_filter] - filter_a[n] * x_u_edge[n,index_filter]
        else:
            index_filter = len(x_l_edge[n,:])-1
            filter_a[n] = (y_l_edge[n,index_filter] - y_l_edge[n,filter_k[n]]) / (x_l_edge[n,index_filter] - x_l_edge[n,filter_k[n]])
            filter_b[n] = y_l_edge[n,index_filter] - filter_a[n] * x_l_edge[n,index_filter]

    for n in range(3):
        x_u_p = x_u_edge[n,:]
        y_u_p = y_u_edge[n,:]
        for i in range(len(x_u_p)):
           if (x_u_p[i] >= x_hinge):
               x_u_edge[n,i] = (x_u_p[i]-x_hinge) * np.cos(angle) - (y_u_p[i]-y_hinge) * np.sin(angle) + x_hinge
               y_u_edge[n,i] = (x_u_p[i]-x_hinge) * np.sin(angle) + (y_u_p[i]-y_hinge) * np.cos(angle) + y_hinge
        x_l_p = x_l_edge[n,:]
        y_l_p = y_l_edge[n,:]
        for i in range(len(x_l_p)):
            if (x_l_p[i] >= x_hinge):
                x_l_edge[n,i] = (x_l_p[i]-x_hinge) * np.cos(angle) - (y_l_p[i]-y_hinge) * np.sin(angle) + x_hinge
                y_l_edge[n,i] = (x_l_p[i]-x_hinge) * np.sin(angle) + (y_l_p[i]-y_hinge) * np.cos(angle) + y_hinge

    # Smoothing

    for n in range(3):

        if angle >= 0.0:
            points_to_remove = []
            for m in range(filter_k[n]+1,len(x_u_edge[n,:])-1): # can't remove last one
                if (y_u_edge[n,m] < filter_a[n]*x_u_edge[n,m]+filter_b[n]):
                    points_to_remove.append(m)
            f_u = interpolate.interp1d(np.delete(x_u_edge[n,:],points_to_remove,None),np.delete(y_u_edge[n,:],points_to_remove,None))
            f_l = interpolate.interp1d(x_l_edge[n,:],y_l_edge[n,:])

        else:
            f_u = interpolate.interp1d(x_u_edge[n,:],y_u_edge[n,:])
            points_to_remove = []
            for m in range(filter_k[n]+1,len(x_l_edge[n,:])-1): # can't remove last one
                if (y_l_edge[n,m] > filter_a[n]*x_l_edge[n,m]+filter_b[n]):
                    points_to_remove.append(m)
            f_l = interpolate.interp1d(np.delete(x_l_edge[n,:],points_to_remove,None),np.delete(y_l_edge[n,:],points_to_remove,None))

        x_distrib = x_canvas-pc_hinge
        for k in range(len(x_distrib)):
            if (x_distrib[k] > 0):
                x_distrib[k] = x_distrib[k]/(1.0-pc_hinge)*(np.max(x_u_edge[n,:])-x_hinge)+x_hinge
            else:
                n_last = k # is the last point before hinge location assuming x_distrib is in increasing order
                x_distrib[k] = x_distrib[k]/pc_hinge*(x_hinge-np.min(x_u_edge[n,:]))+x_hinge
        x_distrib[0] = x_u_edge[n,0]
        x_distrib[len(x_distrib)-1] = x_u_edge[n,len(x_distrib)-1]
        x_u_edge[n,:] = x_distrib
        y_u_edge[n,:] = f_u(x_u_edge[n,:])
        x_l_edge[n,:] = x_distrib
        y_l_edge[n,:] = f_l(x_l_edge[n,:])

    # Redistribute points

    n_points = 65
    n_first = 15

    f_u = interpolate.interp1d(x_u_root,y_u_root)
    x_u_root = np.append(np.append(x_u_root[0:n_first],np.linspace(x_u_root[n_first],x_u_root[-n_last],n_points-n_first-n_last+1)),x_u_root[-n_last+1:len(x_u_root)])
    y_u_root = f_u(x_u_root)
    f_l = interpolate.interp1d(x_l_root,y_l_root)
    x_l_root = np.append(np.append(x_l_root[0:n_first],np.linspace(x_l_root[n_first],x_l_root[-n_last],n_points-n_first-n_last+1)),x_l_root[-n_last+1:len(x_l_root)])
    y_l_root = f_l(x_l_root)

    f_u = interpolate.interp1d(x_u_in_cut,y_u_in_cut)
    x_u_in_cut = np.append(np.append(x_u_in_cut[0:n_first],np.linspace(x_u_in_cut[n_first],x_u_in_cut[-n_last],n_points-n_first-n_last+1)),x_u_in_cut[-n_last+1:len(x_u_in_cut)])
    y_u_in_cut = f_u(x_u_in_cut)
    f_l = interpolate.interp1d(x_l_in_cut,y_l_in_cut)
    x_l_in_cut = np.append(np.append(x_l_in_cut[0:n_first],np.linspace(x_l_in_cut[n_first],x_l_in_cut[-n_last],n_points-n_first-n_last+1)),x_l_in_cut[-n_last+1:len(x_l_in_cut)])
    y_l_in_cut = f_l(x_l_in_cut)


    x_u_edge_new = np.zeros((3,n_points))
    y_u_edge_new = np.zeros((3,n_points))
    x_l_edge_new = np.zeros((3,n_points))
    y_l_edge_new = np.zeros((3,n_points))

    for n in range(3):
        f_u = interpolate.interp1d(x_u_edge[n,:],y_u_edge[n,:])
        x_u_edge_new[n,:] = np.append(np.append(x_u_edge[n,0:n_first],np.linspace(x_u_edge[n,n_first],x_u_edge[n,-n_last],n_points-n_first-n_last+1)),x_u_edge[n,-n_last+1:len(x_u_edge[n,:])])
        y_u_edge_new[n,:] = f_u(x_u_edge_new[n,:])
        f_l = interpolate.interp1d(x_l_edge[n,:],y_l_edge[n,:])
        x_l_edge_new[n,:] = np.append(np.append(x_l_edge[n,0:n_first],np.linspace(x_l_edge[n,n_first],x_l_edge[n,-n_last],n_points-n_first-n_last+1)),x_l_edge[n,-n_last+1:len(x_l_edge[n,:])])
        y_l_edge_new[n,:] = f_l(x_l_edge_new[n,:])

    x_u_edge = x_u_edge_new
    y_u_edge = y_u_edge_new
    x_l_edge = x_l_edge_new
    y_l_edge = y_l_edge_new

    # Propagating profiles

    n_profiles = n_profiles_0 + n_profiles_1 + n_profiles_2 - 2
    x_u = np.zeros((n_profiles,len(x_u_edge[0,:])))
    y_u = np.zeros((n_profiles,len(y_u_edge[0,:])))
    x_l = np.zeros((n_profiles,len(x_l_edge[0,:])))
    y_l = np.zeros((n_profiles,len(y_l_edge[0,:])))

    for n in range(n_profiles_0):
        coeff = n/float(n_profiles_0-1)
        for i in range(len(x_u_root)):
            x_u[n,i] = x_u_root[i] + coeff*(x_u_in_cut[i]-x_u_root[i])
            y_u[n,i] = y_u_root[i] + coeff*(y_u_in_cut[i]-y_u_root[i])
            x_l[n,i] = x_l_root[i] + coeff*(x_l_in_cut[i]-x_l_root[i])
            y_l[n,i] = y_l_root[i] + coeff*(y_l_in_cut[i]-y_l_root[i])

    for n in range(1,n_profiles_1):
        coeff = n/float(n_profiles_1-1)
        for i in range(len(x_u_edge[0,:])):
            x_u[n_profiles_0-1+n,i] = x_u_edge[0,i] + coeff*(x_u_edge[1,i]-x_u_edge[0,i])
            y_u[n_profiles_0-1+n,i] = y_u_edge[0,i] + coeff*(y_u_edge[1,i]-y_u_edge[0,i])
            x_l[n_profiles_0-1+n,i] = x_l_edge[0,i] + coeff*(x_l_edge[1,i]-x_l_edge[0,i])
            y_l[n_profiles_0-1+n,i] = y_l_edge[0,i] + coeff*(y_l_edge[1,i]-y_l_edge[0,i])

    for n in range(1,n_profiles_2):
        coeff = n/float(n_profiles_2-1)
        for i in range(len(x_u_edge[1,:])):
            x_u[n_profiles_0+n_profiles_1-2+n,i] = x_u_edge[1,i] + coeff*(x_u_edge[2,i]-x_u_edge[1,i])
            y_u[n_profiles_0+n_profiles_1-2+n,i] = y_u_edge[1,i] + coeff*(y_u_edge[2,i]-y_u_edge[1,i])
            x_l[n_profiles_0+n_profiles_1-2+n,i] = x_l_edge[1,i] + coeff*(x_l_edge[2,i]-x_l_edge[1,i])
            y_l[n_profiles_0+n_profiles_1-2+n,i] = y_l_edge[1,i] + coeff*(y_l_edge[2,i]-y_l_edge[1,i])



    # Writing files

    for n in range(n_profiles):
        out_file = 'profile_%s_%d.dat' % (case,(n+1))
        out = open(out_file,"w")
        for k in range(len(x_u[n,:])):
            out.write(str(x_u[n,len(x_u[n,:])-1-k]) + ' ' + str(y_u[n,len(x_u[n,:])-1-k]) + '\n')
        for k in range(1,len(x_l[n,:])):
            out.write(str(x_l[n,k]) + ' ' + str(y_l[n,k]) + '\n')
        out.close()

    shutil.copy(SPACE_RUN + '/SPACE/util/profiles/profile_cs.dat', 'profile_cs.dat')

    return wing_width_section_1, wing_width_section_2


