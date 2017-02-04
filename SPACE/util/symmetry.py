#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def mesh2tri(config):

    # Input

    nodes_baseline, elems_baseline = read_mesh(config.FLUID_SURFACE + '.mesh')
    nodes_farfield, elems_farfield = read_mesh(config.FARFIELD_FILENAME)

    # Process

    edges_baseline, sym_nodes_baseline = process(nodes_baseline, elems_baseline)
    edges_farfield, sym_nodes_farfield = process(nodes_farfield, elems_farfield)

    # Output

    out = open(config.FLUID_SURFACE + '.poly',"w")

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

def tri2mesh(config):

    nodes = open(config.FLUID_SURFACE + '.1.node')
    elems = open(config.FLUID_SURFACE + '.1.ele')
    out = open(config.SYMMETRY_FILENAME,"w")

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

