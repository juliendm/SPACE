#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, time, sys, shutil, copy
import numpy as np

from .. import io  as spaceio

# ----------------------------------------------------------------------
#  Load Simulation
# ----------------------------------------------------------------------

def load(config, loadFactor, gravity_vector, safetyFactor_thrust = 1.0, safetyFactor_inertial = 1.0, safetyFactor_non_inertial = 1.0):

    # local copy
    konfig = copy.deepcopy(config)


    thrust_vector = np.array([-float(konfig.THRUST), 0.0, 0.0])*safetyFactor_thrust
    thrust_balance = 0.5;

    nDim = 3
    nNode = 4

    # Read Files
    coord_bdf, elem_bdf, elem_tag_bdf, descriptions = read_bdf(konfig, nDim, nNode)  # OML + Structure
    coord, elem, bdf_corresp                        = read_mesh(konfig, nDim, nNode) # OML only
    pressureCoeff, frictionCoeff                    = read_sol(konfig, nDim, nNode)  # OML only

    # Transpose to BDF
    pressureCoeff_bdf = [0.0 for iPoint_bdf in range(len(coord_bdf))]
    frictionCoeff_bdf = [[0.0 for iDim in range(nDim)] for iPoint_bdf in range(len(coord_bdf))]
    for iPoint in range(len(coord)):
        pressureCoeff_bdf[bdf_corresp[iPoint]] = pressureCoeff[iPoint]
        for iDim in range(nDim):
            frictionCoeff_bdf[bdf_corresp[iPoint]][iDim] = frictionCoeff[iPoint][iDim]

    # Compute Area and Normal
    area, normal = areas_normals(nDim, nNode, coord, elem)
    area_bdf = [0.0 for iPoint_bdf in range(len(coord_bdf))]
    normal_bdf = [[0.0 for iDim in range(nDim)] for iPoint_bdf in range(len(coord_bdf))]
    for iPoint in range(len(coord)):
        area_bdf[bdf_corresp[iPoint]] = area[iPoint]
        for iDim in range(nDim):
            normal_bdf[bdf_corresp[iPoint]][iDim] = normal[iPoint][iDim]

    # sum_x = 0.0
    # sum_y = 0.0
    # sum_z = 0.0
    # sum_area = 0.0
    # for iPoint_bdf in range(nPoint_bdf):
    #     sum_x += normal_bdf[iPoint_bdf][0]
    #     sum_y += normal_bdf[iPoint_bdf][1]
    #     sum_z += normal_bdf[iPoint_bdf][2]
    #     sum_area += area_bdf[iPoint_bdf]
    # print sum_x
    # print sum_y
    # print sum_z
    # print sum_area

    # Compute Load





    # # Control Forces

    # aero_forces = [0.0 for iDim in range(nDim)]
    # for iPoint_bdf in range(nPoint_bdf):

    #     # Pressure
        
    #     pressure = float(konfig.P_DYN_INF)*pressureCoeff_bdf[iPoint_bdf] # + float(konfig.P_INF) # NOT ADDING p_inf CAUSE SPACEPLANE IS NOT PRESSURIZED
    #     for iDim in range(nDim):
    #         aero_forces[iDim] += normal_bdf[iPoint_bdf][iDim]*pressure

    #     # Friction
        
    #     for iDim in range(nDim):
    #         shear_stress = float(konfig.P_DYN_INF)*frictionCoeff_bdf[iPoint_bdf][iDim]
    #         aero_forces[iDim] += area_bdf[iPoint_bdf]*shear_stress

    # for iDim in range(nDim):
    #     aero_forces[iDim] /= float(konfig.P_DYN_INF) * float(konfig.REF_AREA)

    # print "AERO FORCES: ", aero_forces

    # aoa = float(konfig.AoA)*np.pi/180.0

    # lift_coeff = ( -aero_forces[0]*np.sin(aoa) + aero_forces[1]*np.cos(aoa) )
    # drag_coeff = ( aero_forces[0]*np.cos(aoa) + aero_forces[1]*np.sin(aoa) )

    # print lift_coeff, drag_coeff



    nPoint_bdf = len(coord_bdf)
    load_bdf = [[0.0 for iDim in range(nDim)] for iPoint_bdf in range(nPoint_bdf)]

    # NON INERTIAL

    apply_thrust, apply_fuse_r = get_apply_noninertial(descriptions, elem_bdf, elem_tag_bdf)

    for iPoint_bdf in range(nPoint_bdf):

        # AERO

        if not iPoint_bdf in apply_fuse_r:

            # Pressure
            
            pressure = float(konfig.P_DYN_INF)*pressureCoeff_bdf[iPoint_bdf] # + float(konfig.P_INF) # NOT ADDING p_inf CAUSE SPACEPLANE IS NOT PRESSURIZED
            for iDim in range(nDim):
                load_bdf[iPoint_bdf][iDim] += normal_bdf[iPoint_bdf][iDim]*pressure * safetyFactor_non_inertial

            # Friction
            
            for iDim in range(nDim):
                shear_stress = float(konfig.P_DYN_INF)*frictionCoeff_bdf[iPoint_bdf][iDim]
                load_bdf[iPoint_bdf][iDim] += area_bdf[iPoint_bdf]*shear_stress * safetyFactor_non_inertial

        # THRUST

        for iDim in range(nDim):

            if iPoint_bdf in apply_thrust:
                index_thrust = apply_thrust.index(iPoint_bdf)
                if index_thrust == 0:
                    load_bdf[iPoint_bdf][iDim] += thrust_balance*thrust_vector[iDim]
                elif index_thrust == 1:
                    load_bdf[iPoint_bdf][iDim] += (1.0-thrust_balance)*thrust_vector[iDim]


    # INERTIAL

    additional_mass_bdf = get_additional_mass(konfig, descriptions, coord_bdf, elem_bdf, elem_tag_bdf, area_bdf, 15e3, 1e6, 1.0, 1e6)

    for iPoint_bdf in range(nPoint_bdf):
        for iDim in range(nDim):
            load_bdf[iPoint_bdf][iDim] += additional_mass_bdf[iPoint_bdf]*gravity_vector[iDim]*loadFactor*safetyFactor_inertial


    # Write load
    write_load(konfig.LOAD_FILENAME,load_bdf,coord_bdf,elem_bdf)

    # Write Check
    write_check(load_bdf,coord_bdf,elem_bdf)

    # info out
    info = spaceio.State()
    info.FILES.LOAD = konfig.LOAD_FILENAME
    return info

#: def load()







def write_load(filename,load_bdf,coord_bdf,elem_bdf):

    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    load = open(filename,'w')
    load.write(str(nPoint_bdf) + " " + str(nElem_bdf) + "\n")
    for iPoint_bdf in range(nPoint_bdf):
        load.write(str(coord_bdf[iPoint_bdf][0]) + " " + str(coord_bdf[iPoint_bdf][1]) + " " + str(coord_bdf[iPoint_bdf][2]) + " " + str(load_bdf[iPoint_bdf][0]) + " " + str(load_bdf[iPoint_bdf][1]) + " " + str(load_bdf[iPoint_bdf][2]) + "\n")
    for iElem_bdf in range(nElem_bdf):
        load.write(str(elem_bdf[iElem_bdf][0]-1) + " " + str(elem_bdf[iElem_bdf][1]-1)  + " " + str(elem_bdf[iElem_bdf][2]-1) + " " + str(elem_bdf[iElem_bdf][3]-1) + "\n")
    load.close()

#: def write_load()

def write_check(load_bdf,coord_bdf,elem_bdf):

    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    load_mesh = open('struct_load.mesh', 'w')
    load_mesh.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n' + str(nPoint_bdf) + '\n\n')
    for iPoint_bdf in range(nPoint_bdf):
        load_mesh.write(str(coord_bdf[iPoint_bdf][0]) + " " + str(coord_bdf[iPoint_bdf][1]) + " " + str(coord_bdf[iPoint_bdf][2]) + " " + str(iPoint_bdf+1) + "\n")
    load_mesh.write('\nQuadrilaterals\n' + str(nElem_bdf) + '\n\n')
    for iElem_bdf in range(nElem_bdf):
        load_mesh.write(str(elem_bdf[iElem_bdf][0]) + " " + str(elem_bdf[iElem_bdf][1])  + " " + str(elem_bdf[iElem_bdf][2]) + " " + str(elem_bdf[iElem_bdf][3]) + " 0\n")
    load_mesh.write('\nEnd\n')
    load_mesh.close()

    load_sol = open('struct_load.sol', 'w')
    load_sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(nPoint_bdf) + '\n3 1 1 1\n')
    for iPoint_bdf in range(nPoint_bdf):
        load_sol.write(str(load_bdf[iPoint_bdf][0]) + " " + str(load_bdf[iPoint_bdf][1]) + " " + str(load_bdf[iPoint_bdf][2]) + "\n")
    load_sol.write('\nEnd\n')
    load_sol.close()

#: def write_check()







def read_bdf(config, nDim, nNode):

    # Read bdf
    bdf = open(config.STRUCT + '.bdf')
    coord_bdf = []
    elem_bdf = []
    elem_tag_bdf = []
    descriptions = {}
    for line in bdf:
        data = line.split()
        if (line[0]=="$" and len(data) == 3):
            descriptions[data[2].strip().split('/')[0].upper()] = int(data[1])
        elif (line[0]=="G" and len(data) == 6):
            vec = [float(data[3]), float(data[4].strip('*'))]
        elif (line[0]=="*" and len(data) == 5):
            vec.append(float(data[2]))
            coord_bdf.append(vec)
        elif (line[0]=="C" and len(data) == 7):
            elem_bdf.append([int(data[3]), int(data[4]), int(data[5]), int(data[6])])
            elem_tag_bdf.append(int(data[2]))
    bdf.close()

    return coord_bdf, elem_bdf, elem_tag_bdf, descriptions

#: def read_bdf()

def read_mesh(config, nDim, nNode):

    # Read mesh

    mesh = open(config.STRUCT + '_surface.mesh')
    line = mesh.readline()
    while not line.strip() == 'Vertices':
        line = mesh.readline()
    nPoint = int(mesh.readline())
    mesh.readline()
    coord = [[0.0 for iDim in range(nDim)] for iPoint in range(nPoint)]
    bdf_corresp = [0]*nPoint
    for iPoint in range(nPoint):
        data = mesh.readline().split()
        # DO THE ROTATION + MIRROR
        coord[iPoint][0] = float(data[0])
        coord[iPoint][1] = float(data[2])
        coord[iPoint][2] = float(data[1])
        bdf_corresp[iPoint] = int(data[3])-1
    line = mesh.readline()
    while not isInt(line):
        line = mesh.readline()
    nElem = int(line)
    mesh.readline()
    elem = [[0 for iNode in range(nNode)] for iElem in range(nElem)]
    for iElem in range(nElem):
        data = mesh.readline().split()
        for iNode in range(nNode):
            elem[iElem][iNode] = int(data[iNode])-1
    mesh.close()

    return coord, elem, bdf_corresp

#: def read_mesh()

def read_sol(config, nDim, nNode):

    # Read sol
    sol = open(config.STRUCT + '_surface.sol')
    line = sol.readline()
    while not line.strip() == 'SolAtVertices':
        line = sol.readline()
    nPoint = int(sol.readline())
    sol.readline()
    sol.readline()

    pressureCoeff = [0.0 for iPoint in range(nPoint)]
    frictionCoeff = [[0.0 for iDim in range(nDim)] for iPoint in range(nPoint)]

    for iPoint in range(nPoint):
        data = sol.readline().split()
        pressureCoeff[iPoint] = float(data[0])
        # DO THE ROTATION + MIRROR
        frictionCoeff[iPoint][0] = float(data[0+1])
        frictionCoeff[iPoint][1] = float(data[2+1])
        frictionCoeff[iPoint][2] = float(data[1+1])
    sol.close()

    return pressureCoeff, frictionCoeff

#: def read_sol()


def areas_normals(nDim, nNode, coord, elem):

    nPoint = len(coord)
    nElem = len(elem)

    # Compute Normals

    normal = [[0.0 for iDim in range(nDim)] for iPoint in range(nPoint)]
    area = [0.0 for iPoint in range(nPoint)]

    coordElemCG = [[0.0 for iDim in range(nDim)] for iElem in range(nElem)]
    coordEdgeCG = [0.0 for iDim in range(nDim)]
    vec_a = [0.0 for iDim in range(nDim)] 
    vec_b = [0.0 for iDim in range(nDim)] 

    # coordElemCG

    for iElem in range(nElem):
        for iNode in range(nNode):
            iPoint = elem[iElem][iNode]
            for iDim in range(nDim):
                coordElemCG[iElem][iDim] += coord[iPoint][iDim];

    for iElem in range(nElem):
        for iDim in range(nDim):
            coordElemCG[iElem][iDim] /= nNode * 1.0

    # Normal and Area Voronoi

    nNeighnour = 2
    for iElem in range(nElem):
        for iNode in range(nNode):
            iPoint = elem[iElem][iNode]
            for iNeighbour in range(nNeighnour):
                jNode = (iNode + 1 - nNeighnour * iNeighbour) % nNode
                jPoint = elem[iElem][jNode]
                for iDim in range(nDim):
                    coordEdgeCG[iDim] = 0.5 * (coord[iPoint][iDim] + coord[jPoint][iDim])
                if (iNeighbour == 0):
                    for iDim in range(nDim):
                        vec_a[iDim] = coord[iPoint][iDim]-coordElemCG[iElem][iDim]
                        vec_b[iDim] = coordEdgeCG[iDim]-coordElemCG[iElem][iDim]
                else:
                    for iDim in range(nDim):
                        vec_a[iDim] = coord[iPoint][iDim]-coordEdgeCG[iDim]
                        vec_b[iDim] = coordElemCG[iElem][iDim]-coordEdgeCG[iDim]
                normal[iPoint][0] += 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
                normal[iPoint][1] += -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])
                normal[iPoint][2] += 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])

    for iPoint in range(nPoint):
        area[iPoint] = (normal[iPoint][0]**2.0+normal[iPoint][1]**2.0+normal[iPoint][2]**2.0)**0.5

    return area, normal

#: def areas_normals()


def get_n_members():

    nFrame = 13
    nLongeron = 4

    return nFrame, nLongeron

#: def get_n_members()


def get_apply_noninertial(descriptions, elem_bdf, elem_tag_bdf):

    nFrame, nLongeron = get_n_members()

    # THRUST

    thrust_frames = [10]

    tag_thrust_frames = []
    for i in range(nLongeron-1):
        for j in thrust_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_thrust_frames.append(descriptions[desc])

    tag_longerons = []
    for i in [0]:
        for j in range(nFrame-1):
            for desc in ['MLONG:%02d:3:%02d' % (i,j), 'MLONG:%02d:4:%02d' % (i,j)]:
                if desc in descriptions.keys():
                    tag_longerons.append(descriptions[desc])

    tag_members = []
    tag_skin = []
    for desc in descriptions.keys():
        if desc.split(":")[0] in ['MRIBF','MRIBV','MRIBW','MSPARF','MSPARV','MSPARC','MSPARW','MSTRINGC','MSTRINGW','MSKINC','MFRAME','MLONG']:
            tag_members.append(descriptions[desc])
        else:
            tag_skin.append(descriptions[desc])


    apply_thrust_frames = tag_to_apply(tag_thrust_frames, elem_bdf, elem_tag_bdf)
    apply_longerons     = tag_to_apply(tag_longerons    , elem_bdf, elem_tag_bdf)
    apply_skin          = tag_to_apply(tag_skin         , elem_bdf, elem_tag_bdf)

    apply_thrust = intersection(intersection(apply_thrust_frames, apply_longerons), apply_skin)

    # EXHAUST

    tag_fuse_r = []
    for desc in descriptions.keys():
        if 'FUSE_R' in desc:
            tag_fuse_r.append(descriptions[desc])

    apply_fuse_r = tag_to_apply(tag_fuse_r, elem_bdf, elem_tag_bdf)

    return apply_thrust, apply_fuse_r

#: def get_apply_noninertial()


def get_apply_inertial(descriptions, elem_bdf, elem_tag_bdf):

    nFrame, nLongeron = get_n_members()

    avionics_frames = [0,1]
    lox_frames      = [7,8,9,10]
    kero_frames     = [1,2,3]
    engine_frames   = [10,11,16]
    payload_frames  = [3,4,5,6,7]

    # TAGS

    tag_longerons = []
    for i in [0]:
        for j in range(nFrame-1):
            for desc in ['MLONG:%02d:3:%02d' % (i,j), 'MLONG:%02d:4:%02d' % (i,j)]:
                if desc in descriptions.keys():
                    tag_longerons.append(descriptions[desc])

    tag_longerons_payload = []
    for i in [2]:
        for j in range(nFrame-1):
            for desc in ['MLONG:%02d:3:%02d' % (i,j)]:
                if desc in descriptions.keys():
                    tag_longerons_payload.append(descriptions[desc])

    tag_avionics_frames = []
    tag_lox_frames = []
    tag_kero_frames = []
    tag_engine_frames = []
    tag_payload_frames = []
    for i in range(nLongeron-1):
        for j in avionics_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_avionics_frames.append(descriptions[desc])
        for j in lox_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_lox_frames.append(descriptions[desc])
        for j in kero_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_kero_frames.append(descriptions[desc])
        for j in engine_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_engine_frames.append(descriptions[desc])
        for j in payload_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_payload_frames.append(descriptions[desc])

    tag_wing_body = []
    SKIN_FUSE_L = ['FUSE:BOT','FLAP:LOW','FUSE_F']
    SKIN_WING_L = ['LWING:LOW','LWING_T::1']
    SKINS_L = SKIN_FUSE_L + SKIN_WING_L 
    for desc in descriptions.keys():
        for val in SKINS_L:
            if val in desc:
                tag_wing_body.append(descriptions[desc])


    # APPLY

    apply_wing_body     = tag_to_apply(tag_wing_body, elem_bdf, elem_tag_bdf)
    apply_elevon        = [0]
    apply_body_flap     = [0]
    apply_vertical_tail = [0]


    apply_tps       = join(join(apply_wing_body,apply_elevon),apply_body_flap)
    apply_equip     = apply_wing_body

    apply_longerons         = tag_to_apply(tag_longerons        , elem_bdf, elem_tag_bdf)
    apply_longerons_payload = tag_to_apply(tag_longerons_payload, elem_bdf, elem_tag_bdf)

    apply_avionics_frames = tag_to_apply(tag_avionics_frames, elem_bdf, elem_tag_bdf)

    apply_lox_frames      = tag_to_apply(tag_lox_frames     , elem_bdf, elem_tag_bdf)
    apply_kero_frames     = tag_to_apply(tag_kero_frames    , elem_bdf, elem_tag_bdf)
    apply_engine_frames   = tag_to_apply(tag_engine_frames  , elem_bdf, elem_tag_bdf)

    apply_payload_frames  = tag_to_apply(tag_payload_frames , elem_bdf, elem_tag_bdf)


    apply_gear      = [0]
    apply_hydraulic = [0]
    apply_avionics  = intersection(apply_avionics_frames, apply_longerons)
    apply_elec      = apply_avionics

    apply_lox       = intersection(apply_lox_frames     , apply_longerons)
    apply_kero      = intersection(apply_kero_frames    , apply_longerons)
    apply_engine    = intersection(apply_engine_frames  , apply_longerons)

    apply_payload   = intersection(apply_payload_frames , apply_longerons_payload)

    return apply_gear, apply_hydraulic, apply_avionics, apply_elec, apply_equip, apply_lox, apply_kero, apply_engine, apply_tps, apply_payload, apply_wing_body, apply_elevon, apply_body_flap, apply_vertical_tail

#: def get_apply_inertial()

def tag_to_apply(tag, elem_bdf, elem_tag_bdf):

    apply_list = []

    for iElem_bdf in range(len(elem_bdf)):
        if elem_tag_bdf[iElem_bdf] in tag:
            apply_list += (np.array(elem_bdf[iElem_bdf])-1).tolist()

    return np.unique(apply_list).tolist()

#: def tag_to_apply()

def intersection(vec_a, vec_b):

    return [val for val in vec_a if val in vec_b]

#: def intersection()

def join(vec_a, vec_b):

    return np.unique(vec_a + vec_b).tolist()

#: def join()

def get_vehicle_dimensions(coord):

    coord_array = np.array(coord)

    body_length = 18.0 # abs(max(coord_array[:,0])-min(coord_array[:,0])) # FIXED TO NOT INCLUDE VARIATION DUE TO BODYFLAP !!!!!!!!!!
    wing_span = 2.0*abs(max(coord_array[:,2]))

    return body_length, wing_span

#: get_vehicle_dimesions()

def get_apply_surface(apply_list, area_bdf):

    surface = 0.0

    for iPoint_bdf in apply_list:
        surface += area_bdf[iPoint_bdf]

    return surface

#: get_apply_surface()


def get_additional_mass(config, descriptions, coord_bdf, elem_bdf, elem_tag_bdf, area_bdf, half_mass_kg_current_step, half_thrust_newtons, fuel_percentage, pdyn_max):

    nFrame, nLongeron = get_n_members()

    apply_gear, apply_hydraulic, apply_avionics, apply_elec, apply_equip, apply_lox, apply_kero, apply_engine, apply_tps, apply_payload, apply_wing_body, apply_elevon, apply_body_flap, apply_vertical_tail = get_apply_inertial(descriptions, elem_bdf, elem_tag_bdf)

    mass_gear, mass_hydraulic, mass_avionics, mass_elec, mass_equip, mass_lox, mass_kero, mass_engine = get_half_additional_mass_kg_next_step(half_mass_kg_current_step, half_thrust_newtons, fuel_percentage, pdyn_max, coord_bdf, area_bdf, apply_wing_body, apply_elevon, apply_body_flap, apply_vertical_tail)

    additional_mass_bdf = [0.0 for iPoint_bdf in range(len(coord_bdf))]

    surface_equip = get_apply_surface(apply_equip, area_bdf)

    for iPoint_bdf in range(len(coord_bdf)):

        # Landing Gear

        if iPoint_bdf in apply_gear:

            additional_mass_bdf[iPoint_bdf] += mass_gear/len(apply_gear)

        # Hydraulic

        if iPoint_bdf in apply_hydraulic:

            additional_mass_bdf[iPoint_bdf] += mass_hydraulic/len(apply_hydraulic)

        # Avionics

        if iPoint_bdf in apply_avionics:

            additional_mass_bdf[iPoint_bdf] += mass_avionics/len(apply_avionics)

        # Electrical System

        if iPoint_bdf in apply_elec:

            additional_mass_bdf[iPoint_bdf] += mass_elec/len(apply_elec)

        # Equipment

        if iPoint_bdf in apply_equip:

            additional_mass_bdf[iPoint_bdf] += area_bdf[iPoint_bdf]*mass_equip/surface_equip

        # Tank LOX

        if iPoint_bdf in apply_lox:

            additional_mass_bdf[iPoint_bdf] += mass_lox/len(apply_lox)

        # Tank KERO

        if iPoint_bdf in apply_kero:

            additional_mass_bdf[iPoint_bdf] += mass_kero/len(apply_kero)

        # Engine

        if iPoint_bdf in apply_engine:

            additional_mass_bdf[iPoint_bdf] += mass_engine/len(apply_engine)

        # Thermal Protection System

        if iPoint_bdf in apply_tps:

            # https://science.ksc.nasa.gov/shuttle/technology/sts-newsref/sts-tps.html
            # 9 pounds per cubic foot = 144.166 kg/m^3
            # thickness from 1 inch to 5 inches -> 3 inches = 0.0762 m
            # generally, the HRSI tiles are thicker at the forward areas of the orbiter and thinner toward the aft end

            density_tps = 144.166 # kg/m^3
            thickness_tps = 0.0762 # m
            additional_mass_bdf[iPoint_bdf] += area_bdf[iPoint_bdf]*thickness_tps*density_tps

        # Payload

        if iPoint_bdf in apply_payload:

            additional_mass_bdf[iPoint_bdf] += float(config.PAYLOAD_MASS)/len(apply_payload)

    return additional_mass_bdf

#: def get_additional_mass()





def get_half_additional_mass_kg_next_step(half_mass_kg_current_step, half_thrust_newtons, fuel_percentage, pdyn_max, coord_bdf, area_bdf, apply_wing_body, apply_elevon, apply_body_flap, apply_vertical_tail):

    pounds_to_kg = 0.453592
    newtons_to_pounds = 0.224809
    kg_to_pounds = 2.20462
    meters_to_feet = 3.28084

    weight_pounds_current_step = 2.0*half_mass_kg_current_step*kg_to_pounds
    thrust_pounds = 2.0*half_thrust_newtons*newtons_to_pounds

    body_length, wing_span = get_vehicle_dimensions(coord_bdf)
    body_length_feet = body_length*meters_to_feet
    wing_span_feet = wing_span*meters_to_feet

    wing_body_surface_square_feet     = 2.0*get_apply_surface(apply_wing_body, area_bdf)*meters_to_feet*meters_to_feet
    elevon_surface_square_feet        = 2.0*get_apply_surface(apply_elevon, area_bdf)*meters_to_feet*meters_to_feet
    body_flap_surface_square_feet     = 2.0*get_apply_surface(apply_body_flap, area_bdf)*meters_to_feet*meters_to_feet
    vertical_tail_surface_square_feet = 2.0*get_apply_surface(apply_vertical_tail, area_bdf)*meters_to_feet*meters_to_feet

    N_engines = 1
    rocket_expansion_ratio = 0.0

    rho_tank = 0.0
    V_fuel = 0.0             *2.0
    insulation = 0.0 
    S_tank = 0.0             *2.0


    # Landing Gear Weight

    weight_gear = 0.00916*weight_pounds_current_step**1.124
    half_mass_gear = weight_gear*0.5*pounds_to_kg

    # Hydraulic Weight

    weight_hydraylic = 2.64 * ( ( (wing_body_surface_square_feet + elevon_surface_square_feet + body_flap_surface_square_feet + vertical_tail_surface_square_feet)*pdyn_max/1000.0)**0.334 * (body_length_feet + wing_span_feet)**0.5 )
    half_mass_hydraulic = weight_hydraylic*0.5*pounds_to_kg

    # Avionics Weight

    weight_avionics = 66.37*weight_pounds_current_step**0.361
    half_mass_avionics = weight_avionics*0.5*pounds_to_kg

    # Electrical System Weight

    weight_elec = 1.167*weight_pounds_current_step**0.5*body_length_feet**0.25
    half_mass_elec = weight_elec*0.5*pounds_to_kg

    # Equipment Weight

    weight_equip = 5000.0 + 0.01*weight_pounds_current_step ############## Maybe reduce fixed value
    half_mass_equip = weight_equip*0.5*pounds_to_kg

    # Tank LOX Weight

    weight_tank_lox = rho_tank*V_fuel + insulation*S_tank
    half_mass_lox = weight_tank_lox*0.5*pounds_to_kg          # fuel_percentage

    # Tank KERO Weight

    weight_tank_kero = rho_tank*V_fuel + insulation*S_tank
    half_mass_kero = weight_tank_kero*0.5*pounds_to_kg        # fuel_percentage

    # Engine Weight

    weight_engine = 0.00766*thrust_pounds + 0.00033*thrust_pounds*rocket_expansion_ratio**0.5 + 130.0*N_engines
    half_mass_engine = weight_engine*0.5*pounds_to_kg

    return half_mass_gear, half_mass_hydraulic, half_mass_avionics, half_mass_elec, half_mass_equip, half_mass_lox, half_mass_kero, half_mass_engine

#: def get_half_additional_mass_kg_next_step()




def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
