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

def load(config):

    # local copy
    konfig = copy.deepcopy(config)

    # Read bdf
    bdf = open(konfig.STRUCT + '.bdf')
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
    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    # Read mesh
    nDim = 3
    nNode = 4
    mesh = open(konfig.STRUCT + '_surface.mesh')
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

    # Read sol
    sol = open(konfig.STRUCT + '_surface.sol')
    line = sol.readline()
    while not line.strip() == 'SolAtVertices':
        line = sol.readline()
    nPoint = int(sol.readline())
    sol.readline()
    sol.readline()
    pressureCoeff_bdf = [0.0 for iPoint_bdf in range(nPoint_bdf)]
    for iPoint in range(nPoint):
        data = sol.readline().split()
        pressureCoeff_bdf[bdf_corresp[iPoint]] = float(data[0])
    sol.close()

    # Compute Normals

    normal_bdf = [[0.0 for iDim in range(nDim)] for iPoint_bdf in range(nPoint_bdf)]

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

    # normal

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
            normal_bdf[bdf_corresp[iPoint]][0] += 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
            normal_bdf[bdf_corresp[iPoint]][1] += -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])
            normal_bdf[bdf_corresp[iPoint]][2] += 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])

    # sum_x = 0.0
    # sum_y = 0.0
    # sum_z = 0.0
    # for iPoint_bdf in range(nPoint_bdf):
    #     sum_x += normal_bdf[iPoint_bdf][0]
    #     sum_y += normal_bdf[iPoint_bdf][1]
    #     sum_z += normal_bdf[iPoint_bdf][2]
    # print sum_x
    # print sum_y
    # print sum_z

    # Compute Load

    nx = float(konfig.ACCELERATION_X)
    ny = float(konfig.ACCELERATION_Y)
    nz = float(konfig.ACCELERATION_Z)

    loadFactor = np.sqrt(nx*nx+ny*ny+nz*nz)
    gravity_vector = -9.81 * np.array([-nx/loadFactor,-nz/loadFactor,-ny/loadFactor]) # Change of Frame: to Structure Frame

    thrust_vector = np.array([-float(konfig.THRUST), 0.0, 0.0])
    thrust_balance = 0.5;
    #thrust_vector = np.array([-154834.16375284, -22797.20321268, 0.0])
    #thrust_balance = 0.0771128158

    payload_vector = gravity_vector*loadFactor*float(konfig.PAYLOAD_MASS)
    propu_vector = gravity_vector*loadFactor*float(konfig.PROPU_MASS)
    lox_vector = gravity_vector*loadFactor*float(konfig.LOX_MASS)
    kero_vector = gravity_vector*loadFactor*float(konfig.KERO_MASS)
    gnc_vector = gravity_vector*loadFactor*float(konfig.GNC_MASS)

    nLongeron = 4
    nFrame = 13

    thrust_frames = [10]
    gnc_frames = [0,1]
    kero_frames = [1,2,3]
    payload_frames = [3,4,5,6,7]
    lox_frames = [7,8,9,10]
    propu_frames = [10,11,16]

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

    tag_thrust_frames = []
    tag_gnc_frames = []
    tag_kero_frames = []
    tag_payload_frames = []
    tag_lox_frames = []
    tag_propu_frames = []
    for i in range(nLongeron-1):
        for j in thrust_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_thrust_frames.append(descriptions[desc])
        for j in gnc_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_gnc_frames.append(descriptions[desc])
        for j in kero_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_kero_frames.append(descriptions[desc])
        for j in payload_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_payload_frames.append(descriptions[desc])
        for j in lox_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_lox_frames.append(descriptions[desc])
        for j in propu_frames:
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in descriptions.keys():
                    tag_propu_frames.append(descriptions[desc])

    tag_fuse_r = []
    for desc in descriptions.keys():
        if 'FUSE_R' in desc:
            tag_fuse_r.append(descriptions[desc])

    tag_members = []
    tag_skin = []
    for desc in descriptions.keys():
        if desc.split(":")[0] in ['MRIBF','MRIBV','MRIBW','MSPARF','MSPARV','MSPARC','MSPARW','MSTRINGC','MSTRINGW','MSKINC','MFRAME','MLONG']:
            tag_members.append(descriptions[desc])
        else:
            tag_skin.append(descriptions[desc])

    apply_fuse_r = []
    for iElem_bdf in range(nElem_bdf):
        elem_tag_i = elem_tag_bdf[iElem_bdf]
        elem_i = elem_bdf[iElem_bdf]
        if elem_tag_i in tag_fuse_r:
            apply_fuse_r += elem_i
    apply_fuse_r = np.unique(apply_fuse_r).tolist()

    apply_skin = []
    for iElem_bdf in range(nElem_bdf):
        elem_tag_i = elem_tag_bdf[iElem_bdf]
        elem_i = elem_bdf[iElem_bdf]
        if elem_tag_i in tag_skin:
            apply_skin += elem_i
    apply_skin = np.unique(apply_skin).tolist()

    apply_longerons = []
    apply_longerons_payload = []
    for iElem_bdf in range(nElem_bdf):
        elem_tag_i = elem_tag_bdf[iElem_bdf]
        elem_i = elem_bdf[iElem_bdf]
        if elem_tag_i in tag_longerons:
            apply_longerons += elem_i
        if elem_tag_i in tag_longerons_payload:
            apply_longerons_payload += elem_i
    apply_longerons = np.unique(apply_longerons).tolist()
    apply_longerons_payload = np.unique(apply_longerons_payload).tolist()

    apply_thrust_frames = []
    apply_gnc_frames = []
    apply_kero_frames = []
    apply_payload_frames = []
    apply_lox_frames = []
    apply_propu_frames = []
    for iElem_bdf in range(nElem_bdf):
        elem_tag_i = elem_tag_bdf[iElem_bdf]
        elem_i = elem_bdf[iElem_bdf]
        if elem_tag_i in tag_thrust_frames:
            apply_thrust_frames += elem_i
        if elem_tag_i in tag_gnc_frames:
            apply_gnc_frames += elem_i
        if elem_tag_i in tag_kero_frames:
            apply_kero_frames += elem_i
        if elem_tag_i in tag_payload_frames:
            apply_payload_frames += elem_i
        if elem_tag_i in tag_lox_frames:
            apply_lox_frames += elem_i
        if elem_tag_i in tag_propu_frames:
            apply_propu_frames += elem_i
    apply_thrust_frames = np.unique(apply_thrust_frames).tolist()
    apply_gnc_frames = np.unique(apply_gnc_frames).tolist()
    apply_kero_frames = np.unique(apply_kero_frames).tolist()
    apply_payload_frames = np.unique(apply_payload_frames).tolist()
    apply_lox_frames = np.unique(apply_lox_frames).tolist()
    apply_propu_frames = np.unique(apply_propu_frames).tolist()


    apply_thrust = [val for val in [val for val in apply_thrust_frames if val in apply_longerons] if val in apply_skin]
    apply_gnc = [val for val in apply_gnc_frames if val in apply_longerons]
    apply_kero = [val for val in apply_kero_frames if val in apply_longerons]
    apply_payload = [val for val in apply_payload_frames if val in apply_longerons_payload]
    apply_lox = [val for val in apply_lox_frames if val in apply_longerons]
    apply_propu = [val for val in apply_propu_frames if val in apply_longerons]

    # print apply_fuse_r
    # print apply_thrust
    # print apply_gnc
    # print apply_kero
    # print apply_payload
    # print apply_lox
    # print apply_propu

    # aero_forces = [0.0 for iDim in range(nDim)]
    # for iPoint_bdf in range(nPoint_bdf):
    #     pressure = float(konfig.P_DYN_INF)*pressureCoeff_bdf[iPoint_bdf] + float(konfig.P_INF)
    #     for iDim in range(nDim):
    #         aero_forces[iDim] += normal_bdf[iPoint_bdf][iDim]*pressure
    # print "AERO FORCES: ", aero_forces
    # aoa = 9.959197*np.pi/180.0
    # print  aero_forces[0]*np.cos(aoa) + aero_forces[1]*np.sin(aoa) * 2.0
    # print -aero_forces[0]*np.sin(aoa) + aero_forces[1]*np.cos(aoa) * 2.0

    load_bdf = [[0.0 for iDim in range(nDim)] for iPoint_bdf in range(nPoint_bdf)]

    for iPoint_bdf in range(nPoint_bdf):
        pressure = float(konfig.P_DYN_INF)*pressureCoeff_bdf[iPoint_bdf] # + float(konfig.P_INF) # NOT ADDING p_inf CAUSE SPACEPLANE IS NOT PRESSURIZED
        for iDim in range(nDim):
          if not iPoint_bdf+1 in apply_fuse_r:
              load_bdf[iPoint_bdf][iDim] += normal_bdf[iPoint_bdf][iDim]*pressure
          if iPoint_bdf+1 in apply_thrust:
              index_thrust = apply_thrust.index(iPoint_bdf+1)
              if index_thrust == 0:
                  load_bdf[iPoint_bdf][iDim] += thrust_balance*thrust_vector[iDim]
              elif index_thrust == 1:
                  load_bdf[iPoint_bdf][iDim] += (1.0-thrust_balance)*thrust_vector[iDim]
          if iPoint_bdf+1 in apply_gnc:
              load_bdf[iPoint_bdf][iDim] += gnc_vector[iDim]/len(apply_gnc)
              #print '[', coord_bdf[iPoint_bdf][0], ',', coord_bdf[iPoint_bdf][1], ',', coord_bdf[iPoint_bdf][2], ',',gnc_mass/len(apply_gnc), '],'
          if iPoint_bdf+1 in apply_kero:
              load_bdf[iPoint_bdf][iDim] += kero_vector[iDim]/len(apply_kero)
              #print '[', coord_bdf[iPoint_bdf][0], ',', coord_bdf[iPoint_bdf][1], ',', coord_bdf[iPoint_bdf][2], ',', kero_mass/len(apply_kero), '],'
          if iPoint_bdf+1 in apply_payload:
              load_bdf[iPoint_bdf][iDim] += payload_vector[iDim]/len(apply_payload)
              #print '[', coord_bdf[iPoint_bdf][0], ',', coord_bdf[iPoint_bdf][1], ',', coord_bdf[iPoint_bdf][2], ',', payload_mass/len(apply_payload), '],'
          if iPoint_bdf+1 in apply_lox:
              load_bdf[iPoint_bdf][iDim] += lox_vector[iDim]/len(apply_lox)
              #print '[', coord_bdf[iPoint_bdf][0], ',', coord_bdf[iPoint_bdf][1], ',', coord_bdf[iPoint_bdf][2], ',', lox_mass/len(apply_lox), '],'
          if iPoint_bdf+1 in apply_propu:
              load_bdf[iPoint_bdf][iDim] += propu_vector[iDim]/len(apply_propu)
              #print '[', coord_bdf[iPoint_bdf][0], ',', coord_bdf[iPoint_bdf][1], ',', coord_bdf[iPoint_bdf][2], ',', propu_mass/len(apply_propu), '],'

    # Write load

    load = open(konfig.LOAD_FILENAME,'w')
    load.write(str(nPoint_bdf) + " " + str(nElem_bdf) + "\n")
    for iPoint_bdf in range(nPoint_bdf):
        load.write(str(coord_bdf[iPoint_bdf][0]) + " " + str(coord_bdf[iPoint_bdf][1]) + " " + str(coord_bdf[iPoint_bdf][2]) + " " + str(load_bdf[iPoint_bdf][0]) + " " + str(load_bdf[iPoint_bdf][1]) + " " + str(load_bdf[iPoint_bdf][2]) + "\n")
    for iElem_bdf in range(nElem_bdf):
        load.write(str(elem_bdf[iElem_bdf][0]-1) + " " + str(elem_bdf[iElem_bdf][1]-1)  + " " + str(elem_bdf[iElem_bdf][2]-1) + " " + str(elem_bdf[iElem_bdf][3]-1) + "\n")
    load.close()

    # Write BDF Nastran
    # -----------------

    bdf_nastran = open(konfig.STRUCT + '_nastran.bdf','w')

    bdf = open(konfig.STRUCT + '.bdf')
    for line in bdf:
        if (line.strip() != 'END BULK'):
            bdf_nastran.write(line)
    bdf.close()

    # Material Properties
    write_line(bdf_nastran,'MAT1',r=8)
    write_line(bdf_nastran,'1',r=8)
    write_line(bdf_nastran,konfig.MATERIAL_YOUNG_MODULUS,r=16)
    write_line(bdf_nastran,konfig.MATERIAL_POISSON_RATIO,r=8)
    write_line(bdf_nastran,konfig.MATERIAL_DENSITY,r=16)
    bdf_nastran.write('\n')

    # Shell
    thickness = 1.0
    for desc in descriptions.keys():
        write_line(bdf_nastran,'PSHELL',r=8)
        write_line(bdf_nastran,'%d' % (descriptions[desc]),r=8)
        write_line(bdf_nastran,'1',r=8)
        write_line(bdf_nastran,str(thickness),r=8)
        write_line(bdf_nastran,'1',r=16)
        write_line(bdf_nastran,'1',r=8)
        bdf_nastran.write('\n')

    # Loads
    # -----

    write_line(bdf_nastran,'LOAD',r=8)
    write_line(bdf_nastran,'1',r=8)
    write_line(bdf_nastran,'1.0',r=8)
    write_line(bdf_nastran,'1.0',r=8)
    write_line(bdf_nastran,'2',r=8)
    write_line(bdf_nastran,'1.0',r=8)
    write_line(bdf_nastran,'3',r=8)
    bdf_nastran.write('\n')

    # Force
    for iPoint_bdf in range(nPoint_bdf):
        if (load_bdf[iPoint_bdf][0] != 0.0 and load_bdf[iPoint_bdf][1] != 0.0 and load_bdf[iPoint_bdf][2] != 0.0):
            write_line(bdf_nastran,'FORCE*',r=8)
            write_line(bdf_nastran,'2',r=16)
            write_line(bdf_nastran,'%d' % (iPoint_bdf+1),r=16)
            write_line(bdf_nastran,'0',r=16)
            write_line(bdf_nastran,'1.0',r=16)
            bdf_nastran.write('*F1\n')
            write_line(bdf_nastran,'*F1',r=8)
            write_line(bdf_nastran,'%.8E' % load_bdf[iPoint_bdf][0],r=16)
            write_line(bdf_nastran,'%.8E' % load_bdf[iPoint_bdf][1],r=16)
            write_line(bdf_nastran,'%.8E' % load_bdf[iPoint_bdf][2],r=16)
            bdf_nastran.write('\n')

    # Acceleration
    nx = float(konfig.ACCELERATION_X)
    ny = float(konfig.ACCELERATION_Y)
    nz = float(konfig.ACCELERATION_Z)
    loadFactor = np.sqrt(nx*nx+ny*ny+nz*nz)
    gravityVector = -np.array([-nx,-nz,-ny])/np.sqrt(nx*nx+ny*ny+nz*nz) # Change of Frame: to Structure Frame
    write_line(bdf_nastran,'GRAV',r=8)
    write_line(bdf_nastran,'3',r=8)
    write_line(bdf_nastran,'0',r=8)
    write_line(bdf_nastran,'9.81',r=8)
    write_line(bdf_nastran,'%.3f' % gravityVector[0],r=8)
    write_line(bdf_nastran,'%.3f' % gravityVector[1],r=8)
    write_line(bdf_nastran,'%.3f' % gravityVector[2],r=8)
    bdf_nastran.write('\n')

    # Design Variables
    # ----------------

    min_thickness = 0.001
    max_thickness = 3.0

    for desc in descriptions.keys():
        write_line(bdf_nastran,'DESVAR',r=8)
        write_line(bdf_nastran,'%d' % (descriptions[desc]),r=8)
        write_line(bdf_nastran,'V_' + str(descriptions[desc]),r=8)
        write_line(bdf_nastran,str(thickness),r=8)
        write_line(bdf_nastran,str(min_thickness),r=8)
        write_line(bdf_nastran,str(max_thickness),r=8)
        write_line(bdf_nastran,'1.0',r=8)
        bdf_nastran.write('\n')

    for desc in descriptions.keys():
        write_line(bdf_nastran,'DVPREL1',r=8)
        write_line(bdf_nastran,'%d' % (descriptions[desc]),r=8)
        write_line(bdf_nastran,'PSHELL',r=8)
        write_line(bdf_nastran,'%d' % (descriptions[desc]),r=8)
        write_line(bdf_nastran,'T',r=8)
        # write_line(bdf_nastran,str(min_thickness),r=8)
        # write_line(bdf_nastran,str(max_thickness),r=8)
        # write_line(bdf_nastran,str(0.0),r=8)
        bdf_nastran.write('\n')
        write_line(bdf_nastran,'',r=8)
        write_line(bdf_nastran,'%d' % (descriptions[desc]),r=8)
        write_line(bdf_nastran,str(1.0),r=8)
        bdf_nastran.write('\n')

    # Objective
    # ---------

    dresp_id = 0

    write_line(bdf_nastran,'DRESP1',r=8)
    dresp_id += 1
    write_line(bdf_nastran,str(dresp_id),r=8)
    write_line(bdf_nastran,'W',r=8)
    write_line(bdf_nastran,'WEIGHT',r=8)
    bdf_nastran.write('\n')

    # Constraints
    # -----------

    con_id = 0

    con_id += 1
    write_line(bdf_nastran,'DCONADD',r=8)
    write_line(bdf_nastran,str(con_id),r=8)
    write_line(bdf_nastran,str(con_id+1),r=8)
    # write_line(bdf_nastran,str(con_id+2),r=8)
    bdf_nastran.write('\n')

    # Stresses
    con_id += 1

    var_stress = [9,17]
    for var in var_stress:
        write_line(bdf_nastran,'DRESP1',r=8)
        dresp_id += 1
        write_line(bdf_nastran,str(dresp_id),r=8)
        if (var == 9):
            write_line(bdf_nastran,'VMSTR1',r=8)
        elif (var == 17):
            write_line(bdf_nastran,'VMSTR2',r=8)
        write_line(bdf_nastran,'STRESS',r=8)
        write_line(bdf_nastran,'ELEM',r=16)
        write_line(bdf_nastran,str(var),r=8)
        write_line(bdf_nastran,'',r=8)
        check = 7
        for iElem_bdf in range(nElem_bdf):
            if (check == 8):
                check = 0
                bdf_nastran.write('\n')
                write_line(bdf_nastran,'',r=8)
            check +=1
            write_line(bdf_nastran,str(iElem_bdf+1),r=8)
        if (check != 0):
            bdf_nastran.write('\n')

    write_line(bdf_nastran,'DRESP2',r=8)
    dresp_id += 1
    write_line(bdf_nastran,str(dresp_id),r=8)
    write_line(bdf_nastran,'MAXVM',r=8)
    write_line(bdf_nastran,str(con_id),r=8)
    bdf_nastran.write('\n')
    write_line(bdf_nastran,'',r=8)
    write_line(bdf_nastran,'DRESP1',r=8)
    for iDim in range(len(var_stress)):
        write_line(bdf_nastran,str(dresp_id-len(var_stress)+iDim),r=8)
    bdf_nastran.write('\n')

    write_line(bdf_nastran,'DEQATN',r=8)
    write_line(bdf_nastran,str(con_id),r=8)
    bdf_nastran.write('MAXVM(VMSTR1,VMSTR2)=MAX(VMSTR1,VMSTR2)')
    bdf_nastran.write('\n')

    write_line(bdf_nastran,'DCONSTR',r=8)
    write_line(bdf_nastran,str(con_id),r=8)
    write_line(bdf_nastran,str(dresp_id),r=8)
    write_line(bdf_nastran,'1e-09',r=8)
    write_line(bdf_nastran,'1e7',r=8)
#    write_line(bdf_nastran,konfig.MATERIAL_YIELD_STRENGTH,r=8)
    bdf_nastran.write('\n')

    # Optimization
    # ------------

    write_line(bdf_nastran,'DOPTPRM',r=8) 
    write_line(bdf_nastran,'DESMAX',r=8)  
    write_line(bdf_nastran,'500',r=8)    
    write_line(bdf_nastran,'PENAL',r=8) 
    write_line(bdf_nastran,'0.0',r=8)    
    write_line(bdf_nastran,'CT',r=8) 
    write_line(bdf_nastran,'-.03',r=8)  
    write_line(bdf_nastran,'CTMIN',r=8)    
    write_line(bdf_nastran,'.003',r=8)     
    bdf_nastran.write('\n')
    write_line(bdf_nastran,'',r=8) 
    write_line(bdf_nastran,'CONV1',r=8)  
    write_line(bdf_nastran,'1.-5',r=8)    
    write_line(bdf_nastran,'CONV2',r=8) 
    write_line(bdf_nastran,'1.-20',r=8)    
    write_line(bdf_nastran,'CONVDV',r=8) 
    write_line(bdf_nastran,'1.-6',r=8)  
    write_line(bdf_nastran,'CONVPR',r=8)    
    write_line(bdf_nastran,'1.-5',r=8)
    bdf_nastran.write('\n')

# Default

# CONV1 0.001       Objective relative change 2 iterations
# CONV2 1.0E-20     Objective absolute change 2 iterations
# CONVDV 0.0001     Design variables
# CONVPR 0.001      Properties
# CT -0.03          Constraint tolerance (active if value greater than CT)
# CTMIN 0.003                            (violated if value greater than CTMIN)

# DELB 0.0001       Relative finite difference move parameter
# DELP 0.2          Properties: Fractional change allowed
# DELX 0.5          Design variables: Fractional change allowed

# DESMAX 5

# DPMAX 0.5
# DPMIN 0.01        Properties: Minimum move limit

# DXMAX 1.0
# DXMIN 0.05        Design variables: Minimum move limit

# PENAL 0.0         Improve perfo if starting design is infeasible when e.g. 2.0




    bdf_nastran.write('ENDDATA\n')

    bdf_nastran.close()

    # info out
    info = spaceio.State()
    info.FILES.LOAD = konfig.LOAD_FILENAME
    return info

#: def load()

def write_line(file,line,l=0,r=0):
    if l is not 0:
        n = l - len(line)
        for i in range(n):
            line = ' ' + line
    if r is not 0:
        n = r - len(line)
        for i in range(n):
            line = line + ' '
    file.write(line)

def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
