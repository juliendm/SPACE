#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, time
import numpy as np

from .. import run  as spacerun
from .. import io   as spaceio
from .. import util as spaceutil
from ..eval import model as spacemodel
from ..io   import redirect_folder, redirect_output

# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------

def function(func_name, config, state=None):
    
    # initialize
    state = spaceio.State(state)
    
    # redundancy check
    if not state.FUNCTIONS.has_key(func_name):
        if func_name == 'MISSION': # SPEED_APOGEE
            mission(config, state)
        elif func_name == 'AERODYNAMICS': # DRAG, LIFT
            aerodynamics(config, state)
        elif func_name == 'STRUCTURE': # MASS
            structure(config, state)
        else:
            raise Exception, 'unknown function name, %s' % func_name
    
    # prepare output
    func_out = state.FUNCTIONS
    return copy.deepcopy(func_out)

#: def function()

# ----------------------------------------------------------------------
#  Mission Analysis Function
# ----------------------------------------------------------------------

def mission(config, state=None):

    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = spaceio.State(state)

    # console output
    log_mission = 'log_mission.out'

    # ----------------------------------------------------    
    #  Mission Analysis Solution
    # ----------------------------------------------------    
    
    # redundancy check
    mission_done = all([state.FUNCTIONS.has_key(key) for key in ['SPEED_APOGEE']])

    if not mission_done:

        data = spaceutil.ordered_bunch()
        data.lift_sub = spaceutil.ordered_bunch()
        data.lift_sub.sps = config.DATABASE_LIFT_SUB
        data.lift_sub.ranges = np.array([[0.5,0.95],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
        data.lift_sup = spaceutil.ordered_bunch()
        data.lift_sup.sps = config.DATABASE_LIFT_SUP
        data.lift_sup.ranges = np.array([[1.1,9.0],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
        lift_model = spacemodel.RangedModel(data)

        data = spaceutil.ordered_bunch()
        data.drag_sub = spaceutil.ordered_bunch()
        data.drag_sub.sps = config.DATABASE_DRAG_SUB
        data.drag_sub.ranges = np.array([[0.5,0.95],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
        data.drag_sup = spaceutil.ordered_bunch()
        data.drag_sup.sps = config.DATABASE_DRAG_SUP
        data.drag_sup.ranges = np.array([[1.1,9.0],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
        drag_model = spacemodel.RangedModel(data)

        # files to pull
        files = state.FILES
        pull = []; link = []

        # output redirection
        with redirect_folder('MISSION', pull, link) as push:
            with redirect_output(log_mission):

                # # RUN MISSION ANALYSIS SOLUTION # #
                info = spacerun.mission(config, lift_model, drag_model)
                state.update(info)

    # return output 
    mission = spaceutil.ordered_bunch()
    for key in ['SPEED_APOGEE']:
        if state.FUNCTIONS.has_key(key):
            mission[key] = state.FUNCTIONS[key]

    return mission

#: def mission()

# ----------------------------------------------------------------------
#  Aerodynamics Function
# ----------------------------------------------------------------------

def aerodynamics(config, state=None):
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = spaceio.State(state)
    if not state.FILES.has_key('CONFIG_AERO'):
        state.FILES.CONFIG_AERO = config['CONFIG_AERO_FILENAME']

    # console output
    log_direct = 'log_direct.out'
    
    # ----------------------------------------------------    
    #  Generate Mesh
    # ----------------------------------------------------

    geometry(config,state)

    fluid_mesh(config,state)

    # ----------------------------------------------------    
    #  Direct Solution
    # ----------------------------------------------------    

    # redundancy check
    direct_done = all([state.FUNCTIONS.has_key(key) for key in spaceio.optnames_aero[:9]])
    if not direct_done:
    
        config_aero = spaceio.Config(config.CONFIG_AERO_FILENAME)
        config_aero.MESH_FILENAME = config.FLUID_VOLUME + '.su2'
        config_aero.SURFACE_FLOW_FILENAME = config.FLUID_SURFACE_FLOW
        config_aero.MACH_NUMBER = config.MACH_NUMBER
        config_aero.AoA = config.AoA
        config_aero.SIDESLIP_ANGLE = config.SIDESLIP_ANGLE
        config_aero.REF_ORIGIN_MOMENT_X = config.REF_ORIGIN_MOMENT_X
        config_aero.REF_ORIGIN_MOMENT_Y = config.REF_ORIGIN_MOMENT_Y
        config_aero.REF_ORIGIN_MOMENT_Z = config.REF_ORIGIN_MOMENT_Z
        config_aero.REF_LENGTH_MOMENT = config.REF_LENGTH_MOMENT
        config_aero.REF_AREA = config.REF_AREA

        files = state.FILES
        pull = []; link = []
        
        # files: mesh
        name = files.FLUID_VOLUME_SU2
        link.append(name)

        # output redirection
        with redirect_folder('DIRECT', pull, link) as push:
            with redirect_output(log_direct):     
                
                # # RUN DIRECT SOLUTION # #
                info = spacerun.direct(config_aero)
                #spacerun.restart2solution(config_aero,info)
                state.update(info)
                
                # direct files to push
                name = info.FILES['FLUID_SURFACE_FLOW']
                push.extend([name])

    # return output
    aero = spaceutil.ordered_bunch()
    for key in spaceio.optnames_aero:
        if state.FUNCTIONS.has_key(key):
            aero[key] = state.FUNCTIONS[key]

    return aero

#: def aerodynamics()

# ----------------------------------------------------------------------
#  Structure Function
# ----------------------------------------------------------------------

def structure(config, state=None):

    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = spaceio.State(state)

    # console output
    log_structure = 'log_structure.out'

    # ----------------------------------------------------    
    #  Aerodynamics Solution
    # ----------------------------------------------------    

    aerodynamics(config,state)

    # ----------------------------------------------------    
    #  Structure Solution
    # ----------------------------------------------------    
    
    # redundancy check
    structure_done = all([state.FUNCTIONS.has_key(key) for key in ['MASS']])

    if not structure_done:

        # files to pull
        files = state.FILES
        pull = []; link = []

        # files: mesh
        name = files.FLUID_SURFACE_FLOW
        link.append(name)
        name = files.FLUID_SURFACE_MESH
        link.append(name)
        name = files.STRUCT_SURFACE_MESH
        link.append(name)
        name = files.STRUCT_BDF
        link.append(name)

        # output redirection
        with redirect_folder('STRUCTURE', pull, link) as push:
            with redirect_output(log_structure):
                
                # # RUN STRUCTURE SOLUTION # #
                info = spacerun.structure(config)
                state.update(info)

    # return output 
    struct = spaceutil.ordered_bunch()
    for key in ['MASS']:
        if state.FUNCTIONS.has_key(key):
            struct[key] = state.FUNCTIONS[key]

    return struct

#: def structure()

# ----------------------------------------------------------------------
#  Geometry Function
# ----------------------------------------------------------------------

def geometry(config, state=None):
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = spaceio.State(state)
    
    # console output
    log_geometry = 'log_geometry.out'

    # ----------------------------------------------------    
    #  Geometry Solution
    # ----------------------------------------------------    
    
    # redundancy check
    geometry_done = all( [ state.FILES.has_key(key) for key in ['STRUCT_BDF','STRUCT_MESH','STRUCT_SURFACE_MESH','FLUID_SURFACE_MESH'] ] )

    if not geometry_done:

        # files to pull
        files = state.FILES
        pull = []; link = []

        # output redirection
        with redirect_folder('GEOMETRY', pull, link) as push:
            with redirect_output(log_geometry):
                
                # # RUN GEOMETRY SOLUTION # #
                info = spacerun.geometry(config)
                state.update(info)

                # direct files to push
                name = info.FILES['FLUID_SURFACE_MESH']
                push.extend([name])
                if 'STRUCT_SURFACE_MESH' in info.FILES:
                    push.append(info.FILES['STRUCT_SURFACE_MESH'])
                if 'STRUCT_BDF' in info.FILES:
                    push.append(info.FILES['STRUCT_BDF'])
                if 'STRUCT_MESH' in info.FILES:
                    push.append(info.FILES['STRUCT_MESH'])

    # # return output 
    # geo = spaceutil.ordered_bunch()
    # for key in ['STRUCT_BDF','STRUCT_MESH','STRUCT_SURFACE_MESH','FLUID_SURFACE_MESH']:
    #     if state.FILES.has_key(key):
    #         geo[key] = state.FILES[key]

    # return geo

#: def geometry()

# ----------------------------------------------------------------------
#  Fluid Mesh Function
# ----------------------------------------------------------------------

def fluid_mesh(config, state=None):

    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = spaceio.State(state)
    if not state.FILES.has_key('FLUID_BOUNDARY'):
        state.FILES.FLUID_BOUNDARY = config['FLUID_BOUNDARY_FILENAME']
    if not state.FILES.has_key('CORRESPONDANCE'):
        state.FILES.CORRESPONDANCE = config['CORRESPONDANCE_FILENAME']
    
    # console output
    log_fluid_mesh = 'log_fluid_mesh.out'

    # ----------------------------------------------------    
    #  Fluid Mesh Solution
    # ----------------------------------------------------    
    
    # redundancy check
    fluid_mesh_done = all([state.FILES.has_key(key) for key in ['FLUID_BOUNDARY_UPDATED','FLUID_VOLUME_MESH','FLUID_VOLUME_SU2']])

    if not fluid_mesh_done:

        # files to pull
        files = state.FILES
        pull = []; link = []
        
        # files: mesh
        name = files.FLUID_BOUNDARY
        link.append(name)
        name = files.CORRESPONDANCE
        link.append(name)
        name = files.FLUID_SURFACE_MESH
        link.append(name)
        
        # output redirection
        with redirect_folder('FLUID_MESH', pull, link) as push:
            with redirect_output(log_fluid_mesh):     
                
                # # RUN FLUID MESH SOLUTION # #
                info = spacerun.fluid_mesh(config)
                state.update(info)

                # direct files to push
                name = info.FILES['FLUID_VOLUME_SU2']
                push.extend([name])

    # # return output 
    # fluid = spaceutil.ordered_bunch()
    # for key in ['FLUID_BOUNDARY_UPDATED','FLUID_VOLUME_MESH','FLUID_VOLUME_SU2']:
    #     if state.FILES.has_key(key):
    #         fluid[key] = state.FILES[key]

    # return fluid

#: def fluid_mesh()
