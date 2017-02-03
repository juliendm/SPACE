#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np
from optparse import OptionParser

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")

    (options, args)=parser.parse_args()

    if options.filename == None:
        raise Exception("No config file provided. Use -f flag")
    
    in_file = options.filename
    out_file = 'surface.nas'

    inp = open(in_file)
    out = open(out_file,"w")

    # Read inp Write out
    #---------------------

    # Points    

    line = inp.readline()
    while line[0] != 'V':
        line = inp.readline()

    line = inp.readline()
    nPoint = int(line.split()[0])
 #   line = inp.readline()

    print nPoint

    for i in range(nPoint):
        line = inp.readline()
        data = line.split()
        out.write("GRID  , " + str(i+1) + ",," + data[0] + "," + data[1] + "," + data[2] + "\n")

    # Boundary

    line = inp.readline()
    while line[0] != 'T':
        line = inp.readline()

    line = inp.readline()
    nElem_Bound = int(line.split()[0])
#    line = inp.readline()

    print nElem_Bound

    for i in range(nElem_Bound):
        line = inp.readline()
        data = line.split()
        out.write("CTRIA3,  " + str(i+1) + "," + str(1) + "," + data[0] + "," + data[1] + "," + data[2] + ",,0.000000\n")
 
#: def main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
