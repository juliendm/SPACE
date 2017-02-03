#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, time, sys, shutil, copy
import numpy as np

import operator

from .. import io  as spaceio
from interface import GHS       as SPACE_GHS

# -------------------------------------------------------------------
#  Fluid Mesh Simulation
# -------------------------------------------------------------------

def fluid_mesh ( config ): 

    # local copy
    konfig = copy.deepcopy(config)

    # # Read Correspondance
    # corresp = np.loadtxt(konfig.CORRESPONDANCE_FILENAME, dtype='int')

    # # Read Fluid Surface and Fluid Boundary
    # node_fluid_surface = []
    # fluid_surface = open(konfig.FLUID_SURFACE + '.mesh')
    # line = fluid_surface.readline()
    # while line[0] != 'V':
    #     line = fluid_surface.readline()
    # nPoint = int(fluid_surface.readline().split()[0])
    # fluid_surface.readline()
    # for i in range(nPoint):
    #     data = fluid_surface.readline().split()
    #     node_fluid_surface.append([float(data[0]),float(data[1]),float(data[2])])
    # fluid_surface.close()

    # node_fluid_boundary = []
    # fluid_boundary = open(konfig.FLUID_BOUNDARY_FILENAME)
    # line = fluid_boundary.readline()
    # while line[0] != 'V':
    #     line = fluid_boundary.readline()
    # nPoint = int(fluid_boundary.readline().split()[0])
    # fluid_boundary.readline()
    # for i in range(nPoint):
    #     data = fluid_boundary.readline().split()
    #     node_fluid_boundary.append([float(data[0]),float(data[1]),float(data[2])])
    # remaining = ''
    # line = fluid_boundary.readline()
    # while line:
    #     remaining += line
    #     line = fluid_boundary.readline()
    # fluid_boundary.close()
    # nDim = 3
    # for index in range(len(corresp)):
    #     for iDim in range(nDim):
    #         node_fluid_boundary[corresp[index]-1][iDim] = node_fluid_surface[index][iDim]

    # # Write FLuid Boundary Updated
    # fluid_boundary_updated_filename = konfig.FLUID_BOUNDARY_FILENAME.split('.')
    # fluid_boundary_updated_filename = fluid_boundary_updated_filename[0] + '_updated.' + fluid_boundary_updated_filename[1]
    # fluid_boundary_updated = open(fluid_boundary_updated_filename,'w')
    # nPoint = len(node_fluid_boundary)
    # fluid_boundary_updated.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n' + str(nPoint) + '\n\n')
    # for iPoint in range(nPoint):
    #     fluid_boundary_updated.write("%f %f %f 0\n" % (node_fluid_boundary[iPoint][0], node_fluid_boundary[iPoint][1], node_fluid_boundary[iPoint][2]))
    # fluid_boundary_updated.write(remaining) 
    # fluid_boundary_updated.close()

    # # Run Solution
    # SPACE_GHS(konfig)

    os.system('spider2 -O 1 -in ' + konfig.FLUID_SURFACE + '.mesh -Tri -Ref -from 1:1:10000 -to 1 -out surface_merged.mesh -f -f64 > log_spider.out 2>&1')

    yams = open('surface_merged.yams',"w")
    yams.write('Absolute\nGradation 1.3\nMinSize 0.01\nMaxsize 0.2\nGeomApp 0.001\n')
    yams.close()

    os.system('yams2 -f -O -1 -in surface_merged.mesh -out ' + konfig.FLUID_SURFACE + '.d.mesh > log_yams.out')

    mesh2tri()
    os.system('triangle -p ' + konfig.FLUID_SURFACE + '.poly > log_tri.out')
    tri2mesh()

    os.system('spider2 -O 6 -f -f64 -eps 0.0001 -in ' + konfig.FLUID_SURFACE + '.d.mesh farfield.mesh symmetry.mesh -out boundary.mesh >> log_spider.out 2>&1')

    adap_surf = open('adap.surf',"w")
    adap_surf.write('3    1\n\n')
    adap_surf.close()

    adap_source = open('adap.source',"w")
    adap_source.write('volume boundingbox  -2 20  -1  8  -5 7 h 0.2\n')
    adap_source.close()

    os.system('amg -novol -in boundary.mesh -hgrad 2.0 -out boundary.b.mesh > log_amg.out 2>&1')
    os.system('amg -novol -in boundary.b.mesh -hgrad 1.05 -source adap.source -out boundary.b.1.mesh >> log_amg.out 2>&1')
    os.system('amg -novol -in boundary.b.1.mesh -hgrad 1.05 -source adap.source -out boundary.b.2.mesh >> log_amg.out 2>&1')
    os.system('ghs3d -O 1 -in boundary.b.2.mesh -out volume.meshb > log_ghs.out 2>&1')
    os.system('amg -in volume.meshb -hgrad 1.05 -source adap.source -out ' + konfig.FLUID_VOLUME + '.meshb >> log_amg.out 2>&1')

    os.system('meshutils -O 3 -in ' + konfig.FLUID_VOLUME + '.meshb -out ' + konfig.FLUID_VOLUME + ' > log_meshutil.out')

    # info out
    info = spaceio.State()
    #info.FILES.FLUID_BOUNDARY_UPDATED = fluid_boundary_updated_filename
    info.FILES.FLUID_VOLUME_MESH = konfig.FLUID_VOLUME + '.meshb'
    info.FILES.FLUID_VOLUME_SU2 = konfig.FLUID_VOLUME + '.su2'

    return info

def mesh2tri():

    # Input

    nodes_baseline, elems_baseline = read_mesh('fluid_surface.d.mesh')
    nodes_farfield, elems_farfield = read_mesh('farfield.mesh')

    # Process

    edges_baseline, sym_nodes_baseline = process(nodes_baseline, elems_baseline)
    edges_farfield, sym_nodes_farfield = process(nodes_farfield, elems_farfield)

    # Output

    out = open('fluid_surface.poly',"w")

    index = 0
    out.write('%d %d %d %d\n' % (len(sym_nodes_baseline)+3*len(sym_nodes_farfield), 2, 0, 1))
    corresp_baseline = {}
    for sym_node in sym_nodes_baseline:
        index += 1
        corresp_baseline[sym_node] = index
        out.write('%d %f %f %d\n' % (index, nodes_baseline[sym_node-1][0], nodes_baseline[sym_node-1][2], 2))
    corresp_farfield = {}
    for sym_node in sym_nodes_farfield:
        index += 1
        corresp_farfield[sym_node] = index
        out.write('%d %f %f %d\n' % (index, nodes_farfield[sym_node-1][0], nodes_farfield[sym_node-1][2], 3))
    for sym_node in sym_nodes_farfield:
        index += 1
        out.write('%d %f %f %d\n' % (index, nodes_farfield[sym_node-1][0]/8.0, nodes_farfield[sym_node-1][2]/8.0, 0))
    for sym_node in sym_nodes_farfield:
        index += 1
        out.write('%d %f %f %d\n' % (index, nodes_farfield[sym_node-1][0]/15.0, nodes_farfield[sym_node-1][2]/15.0, 0))

    index = 0
    out.write('%d %d\n' % (len(edges_baseline)+len(edges_farfield), 1))
    for edge in edges_baseline:
        index += 1
        out.write('%d %d %d %d\n' % (index, corresp_baseline[edge[0]], corresp_baseline[edge[1]], 2))
    for edge in edges_farfield:
        index += 1
        out.write('%d %d %d %d\n' % (index, corresp_farfield[edge[0]], corresp_farfield[edge[1]], 3))

    out.write('1\n')
    out.write('1 5 0\n')

    out.close()

 
#: def main()

def read_mesh(filename):

    inp = open(filename)
    line = inp.readline()
    while line[0] != 'V':
        line = inp.readline()
    line = inp.readline()
    nPoint = int(line.split()[0])
    nodes = []
    line = inp.readline()
    while len(line.strip().split()) == 0:
        line = inp.readline()
    for index in range(nPoint):
        data = line.split()
        nodes.append([float(data[0]), float(data[1]), float(data[2])])
        line = inp.readline()
    #line = inp.readline()
    while line[0] != 'T':
        line = inp.readline()
    line = inp.readline()
    nElem = int(line.split()[0])
    elems = []
    line = inp.readline()
    while len(line.strip().split()) == 0:
        line = inp.readline()
    for index in range(nElem):
        data = line.split()
        elems.append([int(data[0]),int(data[1]),int(data[2])])
        line = inp.readline()
    inp.close()

    return nodes, elems

def process(nodes, elems):

    edges = []
    sym_nodes = []
    for elem in elems:
        if abs(nodes[elem[0]-1][1]) < 1e-8 and abs(nodes[elem[1]-1][1]) < 1e-8:
            edges.append([elem[0],elem[1]])
            sym_nodes.append(elem[0])
            sym_nodes.append(elem[1])
        elif abs(nodes[elem[1]-1][1]) < 1e-8 and abs(nodes[elem[2]-1][1]) < 1e-8:
            edges.append([elem[1],elem[2]])
            sym_nodes.append(elem[1])
            sym_nodes.append(elem[2])
        elif abs(nodes[elem[2]-1][1]) < 1e-8 and abs(nodes[elem[0]-1][1]) < 1e-8:
            edges.append([elem[2],elem[0]])
            sym_nodes.append(elem[2])
            sym_nodes.append(elem[0])
    sym_nodes = list(set(sym_nodes))
    sym_nodes.sort()

    return edges, sym_nodes

def tri2mesh():

    nodes = open('fluid_surface.1.node')
    elems = open('fluid_surface.1.ele')
    out = open('symmetry.mesh',"w")

    out.write("\nMeshVersionFormatted \n2\n")
    out.write("\nDimension \n3\n")

    # Points

    out.write("\nVertices\n")

    nPoint = int(nodes.readline().split()[0])
    print nPoint

    out.write(str(nPoint))
    out.write("\n\n")

    for index in range(nPoint):
        data = nodes.readline().split()
        out.write(data[1] + ' 0.0 ' + data[2] + ' 0\n')

    # Elements

    out.write("\nTriangles\n")

    nElem_Bound = int(elems.readline().split()[0])
    print nElem_Bound

    out.write(str(nElem_Bound))
    out.write("\n\n")

    for index in range(nElem_Bound):
        data = elems.readline().split()
        out.write(data[1] + ' ' + data[2] + ' ' + data[3] + ' 3\n')

    out.write("\n\n")
    out.write("End\n")

    nodes.close()
    elems.close()
    out.close()

#: def boundary()
