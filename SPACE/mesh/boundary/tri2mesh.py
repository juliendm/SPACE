#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np
import glob

from optparse import OptionParser

def main():

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


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
