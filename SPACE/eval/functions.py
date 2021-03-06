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

from ..util import DesignVariables
from ..surfpack import Surfpack

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
        elif func_name == 'FLUID_MESH':
            fluid_mesh(config, state)
        elif func_name == 'GEOMETRY':
            geometry(config, state)
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
        data.lift_sub.sps = config.MODEL_LIFT_SUB
        data.lift_sub.ranges = np.array([[0.5,0.95],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
        data.lift_sup = spaceutil.ordered_bunch()
        data.lift_sup.sps = config.MODEL_LIFT_SUP
        data.lift_sup.ranges = np.array([[1.1,9.0],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
        lift_model = spacemodel.RangedModel(data)

        data = spaceutil.ordered_bunch()
        data.drag_sub = spaceutil.ordered_bunch()
        data.drag_sub.sps = config.MODEL_DRAG_SUB
        data.drag_sub.ranges = np.array([[0.5,0.95],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
        data.drag_sup = spaceutil.ordered_bunch()
        data.drag_sup.sps = config.MODEL_DRAG_SUP
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

    kriging = True

    # ----------------------------------------------------    
    #  Generate Geometry
    # ----------------------------------------------------

    geometry(config,state)

    # ----------------------------------------------------    
    #  Compute Aero
    # ----------------------------------------------------

    if kriging:

          return aerodynamics_kriging(config, state)

    else:

        config_aero = spaceio.Config(config.CONFIG_AERO_FILENAME)

        Mach = float(config.MACH_NUMBER)

        if Mach > 1.0:

            config_aero.MGLEVEL= '0'
            config_aero.MGCYCLE= 'W_CYCLE'
            config_aero.EXT_ITER= '3000'

            if Mach < 2.0:                           # 1.0 - 2.0
                config_aero.CFL_NUMBER= '4.0'
            elif Mach < 3.0 :                        # 2.0 - 3.0
                config_aero.CFL_NUMBER= '3.5'
            elif Mach < 5.0:                         # 3.0 - 5.0
                config_aero.CFL_NUMBER= '3.5'
            elif Mach < 8.0:                         # 5.0 - 8.0
                config_aero.CFL_NUMBER= '2.5'
            else:                                    # >= 8.0
                config_aero.CFL_NUMBER= '2.0'

            config_aero.CONV_CRITERIA= 'CAUCHY'
            config_aero.CAUCHY_FUNC_FLOW= 'DRAG'
            config_aero.CAUCHY_EPS= '1E-4'
            config_aero.CAUCHY_ELEMS= '100'


        else:
            
            config_aero.MGLEVEL= '3'
            config_aero.MGCYCLE= 'W_CYCLE'
            config_aero.EXT_ITER= '3000'

            config_aero.CFL_NUMBER= '5.0'

    #        config_aero.MGLEVEL= '3'
    #        config_aero.MGCYCLE= 'W_CYCLE'
    #        config_aero.CFL_NUMBER= '2.0'
    #        config_aero.EXT_ITER= '3000'
    #        config_aero.CFL_ADAPT= 'YES'
    #        config_aero.CFL_ADAPT_PARAM= '( 1.5, 0.5, 2.0, 5.0 )'

            config_aero.CONV_CRITERIA= 'CAUCHY'
            config_aero.CAUCHY_FUNC_FLOW= 'LIFT'
            config_aero.CAUCHY_EPS= '1E-4'
            config_aero.CAUCHY_ELEMS= '100'


        
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

        fluid_mesh(config,state)

        # ----------------------------------------------------    
        #  Direct Solution
        # ----------------------------------------------------  

        # redundancy check
        direct_done = all([state.FUNCTIONS.has_key(key) for key in spaceio.optnames_aero[:9]])

        if not direct_done:
        
            config_aero.MESH_FILENAME = config.FLUID_VOLUME + '.su2'
            config_aero.SURFACE_FLOW_FILENAME = config.FLUID_SURFACE_FLOW

            config_aero.MACH_NUMBER = config.MACH_NUMBER
            config_aero.REYNOLDS_NUMBER = config.REYNOLDS_NUMBER
            config_aero.AoA = config.AoA

            config_aero.SIDESLIP_ANGLE = config.SIDESLIP_ANGLE
            config_aero.REF_ORIGIN_MOMENT_X = config.REF_ORIGIN_MOMENT_X
            config_aero.REF_ORIGIN_MOMENT_Y = config.REF_ORIGIN_MOMENT_Y
            config_aero.REF_ORIGIN_MOMENT_Z = config.REF_ORIGIN_MOMENT_Z
            config_aero.REF_LENGTH_MOMENT = config.REF_LENGTH_MOMENT
            config_aero.REF_AREA = config.REF_AREA

            config_aero.NUMBER_PART = int(config.NUMBER_PART)

            files = state.FILES
            pull = []; link = []
            
            # files: mesh
            name = files.FLUID_VOLUME_MESH
            #name = files.FLUID_VOLUME_SU2
            link.append(name)
            name = files.BOUNDARY_BACK_MESH
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

        else:
            print 'Aerodynamics done'

        # return output
        aero = spaceutil.ordered_bunch()
        for key in spaceio.optnames_aero:
            if state.FUNCTIONS.has_key(key):
                aero[key] = state.FUNCTIONS[key]

        return aero

#: def aerodynamics()

# ----------------------------------------------------------------------
#  Aerodynamics Kriging Function
# ----------------------------------------------------------------------

def aerodynamics_kriging(config, state=None):


    precomputed = False

    # console output
    log_direct = 'log_direct.out'



    # redundancy check

    if not state.FILES.has_key('FLUID_SURFACE_FLOW'):
        fluid_surface_flow_filename = config.FLUID_SURFACE_FLOW + '.sol'
        if os.path.exists(fluid_surface_flow_filename):
            state.FILES.FLUID_SURFACE_FLOW = fluid_surface_flow_filename

    aerodynamics_kriging_done = all( [ state.FILES.has_key(key) for key in ['FLUID_SURFACE_FLOW'] ] )


    if not aerodynamics_kriging_done:

        pull = []; link = []

        # output redirection
        with redirect_folder('DIRECT', pull, link) as push:
            with redirect_output(log_direct):  

                # if Krigings are loaded: aerodynamics_kriging() uses krigings
                if not precomputed:

                    n_models = 2180

                    desvar = DesignVariables()



                    load_models = True

                    if load_models:

                        models_folder = '../../../../RESPONSE_SURFACE_DV_SUP/DESIGNS/MODELS'

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







                    model = Surfpack('placeholder', desvar.ndim)
                    dvs = desvar.pack(config)

                    fluid_surface_flow_filename = config.FLUID_SURFACE_FLOW + '.sol'

                    sol = open(fluid_surface_flow_filename, 'w')
                    sol.write('MeshVersionFormatted 2\n\nDimension 3\n\nSolAtVertices\n' + str(n_models) + '\n4 1 1 1 1\n\n')

                    for index in range(n_models):

                        model.name = 'CP_%05d' % (index+1)
                        cp = model.eval(dvs)

                        model.name = 'CFX_%05d' % (index+1)
                        cfx = model.eval(dvs)

                        model.name = 'CFY_%05d' % (index+1)
                        cfy = model.eval(dvs)

                        model.name = 'CFZ_%05d' % (index+1)
                        cfz = model.eval(dvs)

                        sol.write('%s %s %s %s\n' % (cp, cfx, cfy, cfz))

                    sol.write('\nEnd\n')
                    sol.close()

                    state.FILES.FLUID_SURFACE_FLOW = fluid_surface_flow_filename

                    # direct files to push
                    name = state.FILES.FLUID_SURFACE_FLOW
                    push.extend([name])




                # if Krigings are not loaded: link a file .sol with adequate dsn number
                else:

                    #os.link()

                    raise NotImplementedError

    else:
        print 'Aerodynamics Kriging done'

    # return empty
    return spaceutil.ordered_bunch()

#: def aerodynamics_kriging()

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
    structure_done = False #####################################################################

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

    state.FILES.pop('STRUCT_BDF', None) ####################################

    if not state.FILES.has_key('STRUCT_BDF'):
        struct_bdf_filename = config.STRUCT + '.bdf'
        if os.path.exists(struct_bdf_filename):
            state.FILES.STRUCT_BDF = struct_bdf_filename
    if not state.FILES.has_key('STRUCT_MESH'):
        struct_mesh_filename = config.STRUCT + '.mesh'
        if os.path.exists(struct_mesh_filename):
            state.FILES.STRUCT_MESH = struct_mesh_filename
    if not state.FILES.has_key('STRUCT_SURFACE_MESH'):
        struc_surface_mesh_filename = config.STRUCT + '_surface.mesh'
        if os.path.exists(struc_surface_mesh_filename):
            state.FILES.STRUCT_SURFACE_MESH = struc_surface_mesh_filename
    if not state.FILES.has_key('FLUID_SURFACE_MESH'):
        fluid_surface_mesh_filename = config.FLUID_SURFACE + '.mesh'
        if os.path.exists(fluid_surface_mesh_filename):
            state.FILES.FLUID_SURFACE_MESH = fluid_surface_mesh_filename
    if not state.FILES.has_key('FLUID_SURFACE_BACK_MESH'):
        fluid_surface_back_mesh_filename = config.FLUID_SURFACE + '_back.mesh'
        if os.path.exists(fluid_surface_back_mesh_filename):
            state.FILES.FLUID_SURFACE_BACK_MESH = fluid_surface_back_mesh_filename

    if config.STRUCT == 'NONE':
        geometry_done = all( [ state.FILES.has_key(key) for key in ['FLUID_SURFACE_MESH','FLUID_SURFACE_BACK_MESH'] ] )
    else:
        geometry_done = all( [ state.FILES.has_key(key) for key in ['STRUCT_BDF','STRUCT_MESH','STRUCT_SURFACE_MESH','FLUID_SURFACE_MESH','FLUID_SURFACE_BACK_MESH'] ] )

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
                if 'FLUID_SURFACE_MESH' in info.FILES:
                    # name = info.FILES['FLUID_SURFACE_MESH']
                    # push.extend([name])
                    push.append(info.FILES['FLUID_SURFACE_MESH'])
                if 'FLUID_SURFACE_BACK_MESH' in info.FILES:
                    push.append(info.FILES['FLUID_SURFACE_BACK_MESH'])
                if 'STRUCT_SURFACE_MESH' in info.FILES:
                    push.append(info.FILES['STRUCT_SURFACE_MESH'])
                if 'STRUCT_BDF' in info.FILES:
                    push.append(info.FILES['STRUCT_BDF'])
                if 'STRUCT_MESH' in info.FILES:
                    push.append(info.FILES['STRUCT_MESH'])

    else:
        print 'Geometry done'

    # return output 
    geo = spaceutil.ordered_bunch()
    for key in ['STRUCT_BDF','STRUCT_MESH','STRUCT_SURFACE_MESH','FLUID_SURFACE_MESH','FLUID_SURFACE_BACK_MESH']:
        if state.FILES.has_key(key):
            geo[key] = state.FILES[key]

    return geo

#: def geometry()

# ----------------------------------------------------------------------
#  Fluid Mesh Function
# ----------------------------------------------------------------------

def fluid_mesh(config, state=None):

    config_aero = spaceio.Config(config.CONFIG_AERO_FILENAME)
    config.INITIAL_SPACING = str(initial_spacing(config, config_aero))

    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = spaceio.State(state)

    if not state.FILES.has_key('FARFIELD'):
        state.FILES.FARFIELD = config['FARFIELD_FILENAME']
    
    # console output
    log_fluid_mesh = 'log_fluid_mesh.out'

    # ----------------------------------------------------    
    #  Generate Geometry
    # ----------------------------------------------------

    geometry(config,state)

    # ----------------------------------------------------    
    #  Fluid Mesh Solution
    # ----------------------------------------------------    

    # redundancy check
    fluid_mesh_done = all([state.FILES.has_key(key) for key in ['FLUID_VOLUME_MESH','BOUNDARY_BACK_MESH']])

    if not fluid_mesh_done:

        # files to pull
        files = state.FILES
        pull = []; link = []
        
        # files: mesh
        name = files.FARFIELD
        link.append(name)
        name = files.FLUID_SURFACE_MESH
        link.append(name)
        name = files.FLUID_SURFACE_BACK_MESH
        link.append(name)
        
        # output redirection
        with redirect_folder('FLUID_MESH', pull, link) as push:
            with redirect_output(log_fluid_mesh):     
                
                # # RUN FLUID MESH SOLUTION # #
                info = spacerun.fluid_mesh(config)
                state.update(info)

                # direct files to push
                name = info.FILES['FLUID_VOLUME_MESH']
                push.extend([name])
                push.append(info.FILES['BOUNDARY_BACK_MESH'])



    else:
        print 'Fluid Mesh done'

    # return output 
    fluid = spaceutil.ordered_bunch()
    for key in ['FLUID_VOLUME_MESH','BOUNDARY_BACK_MESH']:
        if state.FILES.has_key(key):
            fluid[key] = state.FILES[key]

    return fluid

#: def fluid_mesh()

def initial_spacing(config, config_aero, yplus=1, BL='turb'):

    Re = float(config.REYNOLDS_NUMBER)
    Ma = float(config.MACH_NUMBER)

    c = (float(config_aero.GAMMA_VALUE) * float(config_aero.GAS_CONSTANT) * float(config_aero.FREESTREAM_TEMPERATURE))**0.5
    mu = float(config_aero.MU_REF) * (float(config_aero.FREESTREAM_TEMPERATURE)/float(config_aero.MU_T_REF))**(3.0/2.0) * (float(config_aero.MU_T_REF) + float(config_aero.SUTHERLAND_CONSTANT))/(float(config_aero.FREESTREAM_TEMPERATURE) + float(config_aero.SUTHERLAND_CONSTANT))
    u = Ma*c
    rho = (Re * mu) / (u * float(config_aero.REYNOLDS_LENGTH))

    if (BL=='turb'):
        Cf = 0.026 / (Re)**(1.0/7.0)
    else:
        Cf = 0.664 / (Re**0.5)

    tau_w = (Cf * rho * u**2) / 2.0
    u_fric = (tau_w / rho)**0.5

    return (yplus * mu) / (u_fric * rho)


