#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np
import glob

from optparse import OptionParser

def main():

    nodes = open('fluid_surface.1.node')
    inp = open('fluid_surface.poly')
    out = open('fluid_surface.a.poly',"w")

    # node

    nNodes = int(nodes.readline().split()[0])
    add = []
    for index in range(nNodes):
        data = nodes.readline().split()
        if (int(data[3]) == 0):
            add.append([float(data[1]), float(data[2])])

    # inp out

    nNodes = int(inp.readline().split()[0])
    out.write('%d 2 0 1\n' % (nNodes + len(add)))
    for index in range(nNodes):
        out.write(inp.readline())
    for index in range(len(add)):
        out.write('%d %f %f 0\n' % (index+1+nNodes, add[index][0], add[index][1]))

    line = inp.readline()
    while line:
        out.write(line)
        line = inp.readline()

    nodes.close()
    inp.close()
    out.close()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
