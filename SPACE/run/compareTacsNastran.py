#!/usr/bin/env python 

import os, time, sys, shutil, copy
import numpy as np
from optparse import OptionParser

import operator

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE import io   as spaceio

sys.path.append(os.environ['STRUCTURE_RUN'])
from mpi4py import MPI
from baseclasses import *
from tacs import *

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main(case):

    if (case == 'NASTRAN_PRE'):

        dvs_file = 'x_final.dat'
        dvs = open(dvs_file)
        x_dvs = np.loadtxt(dvs);
        dvs.close()

        config = spaceio.Config('../config_DSN.cfg')
        konfig = copy.deepcopy(config)

        computeNastran(konfig, x_dvs)

    elif (case == 'NASTRAN_POST'):

        postproNastran('struct_nastran.f06','struct_nastran.sol')

    elif (case == 'TACS'):

        dvs_file = 'x_final.dat'
        dvs = open(dvs_file)
        x_dvs = np.loadtxt(dvs);
        dvs.close()

        config = spaceio.Config('../config_DSN.cfg')
        konfig = copy.deepcopy(config)

        computeTacs(konfig, x_dvs)

#: def main()

def computeTacs(config, x_dvs):

    gcomm = comm = MPI.COMM_WORLD

    # Material properties

    material_rho = float(config.MATERIAL_DENSITY)
    material_E = float(config.MATERIAL_YOUNG_MODULUS)
    material_ys = 2e7 #float(config.MATERIAL_YIELD_STRENGTH)
    material_nu = float(config.MATERIAL_POISSON_RATIO)
    kcorr = 5.0/6.0

    t = 1.0
    tMin = 0.0016 # 0.0016
    tMax = 3.0 # 0.020

    nx = float(config.ACCELERATION_X)
    ny = float(config.ACCELERATION_Y)
    nz = float(config.ACCELERATION_Z)
    loadFactor = np.sqrt(nx*nx+ny*ny+nz*nz)
    gravityVector = -9.81 * np.array([-nx,-nz,-ny])/np.sqrt(nx*nx+ny*ny+nz*nz) # Change of Frame: to Structure Frame

    KSWeight = 80.0
    evalFuncs = ['mass','ks0','mf0']
    SP = StructProblem('lc0', loadFactor=loadFactor, loadFile='../' + config.LOAD_FILENAME, evalFuncs=evalFuncs)
    structOptions = {'transferSize':0.5, 'transferGaussOrder':3}
    FEASolver = pytacs.pyTACS(config.STRUCT + '.bdf', comm=comm, options=structOptions)

    ncoms = FEASolver.nComp
    for i in range(0,ncoms):
        dv_name = 'stru_'+str(i)
        FEASolver.addDVGroup(dv_name, include = i)


    def conCallBack(dvNum, compDescripts, userDescript, specialDVs, **kargs):
        con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, t, dvNum, tMin, tMax)
        scale = [100.0]
        return con, scale
    FEASolver.createTACSAssembler(conCallBack)

    # Mass Functions
    FEASolver.addFunction('mass', functions.StructuralMass)

    # KS Functions
    ks0 = FEASolver.addFunction('ks0', functions.AverageKSFailure, KSWeight=KSWeight, loadFactor=1.0)
    mf0 = FEASolver.addFunction('mf0', functions.AverageMaxFailure)

    # Load Factor
    FEASolver.setOption('gravityVector',gravityVector.tolist())
    FEASolver.addInertialLoad(SP)

    x_dvs_filtered = []
    for x_dv in x_dvs:
        x_dvs_filtered.append(float('%.6f' % x_dv))

    # Write Thickness Sol
    n_point_bdf = 0
    elem_bdf = []
    elem_tag_bdf = []
    bdf = open(config.STRUCT + '.bdf')
    for line in bdf:
        data = line.split()
        if (line[0]=="G" and len(data) == 6):
            n_point_bdf += 1
        elif (line[0]=="C" and len(data) == 7):
            elem_bdf.append([int(data[3]), int(data[4]), int(data[5]), int(data[6])])
            elem_tag_bdf.append(int(data[2]))
    bdf.close()
    n_elem_bdf = len(elem_bdf)
    thickness_point = [0.0 for i_point_bdf in range(n_point_bdf)]
    thickness_point_count = [0 for i_point_bdf in range(n_point_bdf)]
    for i_elem_bdf in range(n_elem_bdf):
        for i_node in range(4):
            thickness_point[elem_bdf[i_elem_bdf][i_node]-1] += x_dvs_filtered[elem_tag_bdf[i_elem_bdf]-1]
            thickness_point_count[elem_bdf[i_elem_bdf][i_node]-1] += 1
    for i_point_bdf in range(n_point_bdf):
        thickness_point[i_point_bdf] = thickness_point[i_point_bdf]/thickness_point_count[i_point_bdf]
    write_sol_1('thicknesses.sol',thickness_point)

    # Compute
    FEASolver.setDesignVars(x_dvs_filtered)
    FEASolver(SP)
    FEASolver.writeMeshDisplacements(SP, "struct_tacs.sol")
    #FEASolver.writeMeshForces(SP, "struct_tacs.sol")
    FEASolver.writeSolution()

def computeNastran(config, x_dvs):

    # Read bdf
    # --------

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
    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    # Write bdf Nastran
    # -----------------

    bdf_nastran = open(config.STRUCT + '_nastran.bdf','w')
    bdf = open(config.STRUCT + '.bdf')
    for line in bdf:
        if (line.strip() != 'END BULK'):
            bdf_nastran.write(line)
    bdf.close()

    # Material Properties
    write_line(bdf_nastran,'MAT1',r=8)
    write_line(bdf_nastran,'1',r=8)
    write_line(bdf_nastran,config.MATERIAL_YOUNG_MODULUS,r=16)
    write_line(bdf_nastran,config.MATERIAL_POISSON_RATIO,r=8)
    write_line(bdf_nastran,config.MATERIAL_DENSITY,r=16)
    bdf_nastran.write('\n')

    # Shell
    for desc in descriptions.keys():
        write_line(bdf_nastran,'PSHELL',r=8)
        write_line(bdf_nastran,'%d' % (descriptions[desc]),r=8)
        write_line(bdf_nastran,'1',r=8)
        #write_line(bdf_nastran,('%.6f' % x_dvs[descriptions[desc]-1])[1:],r=8)
        write_line(bdf_nastran,('%.6f' % 0.0016),r=8)
        write_line(bdf_nastran,'1',r=16)
        write_line(bdf_nastran,'1',r=8)
        bdf_nastran.write('\n')

    # Loads
    # -----

    load = open('../' + config.LOAD_FILENAME)
    nPoint_bdf = int(load.readline().split()[0])
    nDim = 3
    load_bdf = [[0.0 for iDim in range(nDim)] for iPoint_bdf in range(nPoint_bdf)]
    for iPoint_bdf in range(nPoint_bdf):
        data = load.readline().split()
        load_bdf[iPoint_bdf][0] = float(data[3])
        load_bdf[iPoint_bdf][1] = float(data[4])
        load_bdf[iPoint_bdf][2] = float(data[5])
    load.close()

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
    nx = float(config.ACCELERATION_X)
    ny = float(config.ACCELERATION_Y)
    nz = float(config.ACCELERATION_Z)
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
#    write_line(bdf_nastran,config.MATERIAL_YIELD_STRENGTH,r=8)
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

#: def computeNastran()

def postproNastran(input_filename,output_filename):

    coord = []
    elem = []

    f06 = open(input_filename)
    for line in f06:
        data = line.split()
        if (len(data) == 7 and data[1] == 'GRID'):
            index = int(data[2].replace('*',''))
            vec = [float(data[4]),float(data[5])]
        elif (len(data) == 5 and data[1] == '*' and data[3] == '0'):
            vec.append(float(data[2]))
            coord.append(vec)
        elif ((len(data) == 7 or len(data) == 8) and data[1][0] == 'C'):
            index = int(data[2])
            mat = int(data[3])
            vec = [int(data[4]),int(data[5]),int(data[6])]
            if (len(data) == 8): vec.append(int(data[7]))
            elem.append(vec)
    f06.close()

    n_point = len(coord)
    n_elem = len(elem)

    print n_point
    print n_elem

    n_dim = 3
    force = [[0.0 for i_dim in range(n_dim)] for i_point in range(n_point)]
    disp = [[0.0 for i_dim in range(n_dim)] for i_point in range(n_point)]
    rotation = [[0.0 for i_dim in range(n_dim)] for i_point in range(n_point)]

    thickness = [0.0 for i_elem in range(n_elem)]
    vonmises = [[0.0 for index in range(2)] for i_elem in range(n_elem)]

    count = 0
    secondLineVonMises = False
    secondLineForces = False

    f06 = open(input_filename)
    for line in f06:
        data = line.split()
        if (secondLineVonMises and len(data) == 8):
            secondLineVonMises = False
            vonmises[index-1][1] = float(data[7])
        elif (secondLineForces and len(data) == 5):
            secondLineForces = False
            fx = float(data[2]); fy = float(data[3]); fz = float(data[4])
            force[index-1][0] = fx; force[index-1][1] = fy; force[index-1][2] = fz
        elif (len(data) == 7 and data[1] == 'FORCE'):
            index = int(data[3])
            secondLineForces = True
        elif (len(data) == 8 and data[1] == 'G'):
            index = int(data[0])
            dx = float(data[2]); dy = float(data[3]); dz = float(data[4])
            disp[index-1][0] = dx; disp[index-1][1] = dy; disp[index-1][2] = dz
            rx = float(data[5]); ry = float(data[6]); rz = float(data[7])
            rotation[index-1][0] = rx; rotation[index-1][1] = ry; rotation[index-1][2] = rz
        elif (len(data) == 10 and data[0] == '0'):
            index = int(data[1])
            thickness[index-1] = 2.0*abs(float(data[2]))
            vonmises[index-1][0] = float(data[9])
            secondLineVonMises = True
    f06.close()

    #write_sol_3(output_filename,disp)

    thickness_point = [0.0 for i_point in range(n_point)]
    thickness_point_count = [0 for i_point in range(n_point)]
    for i_elem in range(n_elem):
        thickness_point[elem[i_elem][0]-1] += thickness[i_elem]
        thickness_point_count[elem[i_elem][0]-1] += 1
        thickness_point[elem[i_elem][1]-1] += thickness[i_elem]
        thickness_point_count[elem[i_elem][1]-1] += 1
        thickness_point[elem[i_elem][2]-1] += thickness[i_elem]
        thickness_point_count[elem[i_elem][2]-1] += 1
        thickness_point[elem[i_elem][3]-1] += thickness[i_elem]
        thickness_point_count[elem[i_elem][3]-1] += 1
    for i_point in range(n_point):
        thickness_point[i_point] = thickness_point[i_point]/thickness_point_count[i_point]

    vonmises_point = [0.0 for i_point in range(n_point)]
    vonmises_point_count = [0 for i_point in range(n_point)]
    for i_elem in range(n_elem):
        vonmises_point[elem[i_elem][0]-1] += max(vonmises[i_elem][0],vonmises[i_elem][1])
        vonmises_point_count[elem[i_elem][0]-1] += 1
        vonmises_point[elem[i_elem][1]-1] += max(vonmises[i_elem][0],vonmises[i_elem][1])
        vonmises_point_count[elem[i_elem][1]-1] += 1
        vonmises_point[elem[i_elem][2]-1] += max(vonmises[i_elem][0],vonmises[i_elem][1])
        vonmises_point_count[elem[i_elem][2]-1] += 1
        vonmises_point[elem[i_elem][3]-1] += max(vonmises[i_elem][0],vonmises[i_elem][1])
        vonmises_point_count[elem[i_elem][3]-1] += 1
    for i_point in range(n_point):
        vonmises_point[i_point] = vonmises_point[i_point]/vonmises_point_count[i_point]

    write_sol_1(output_filename,thickness_point)

#: def postproNastran()

def write_sol_1(sol_file,solution):
 
    sol = open(sol_file,'w')
    sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(len(solution)) + '\n1 1\n')
    for i_point in range(len(solution)):
        sol.write(str(solution[i_point]) + "\n")
    sol.write('\nEnd\n')
    sol.close()

def write_sol_3(sol_file,solution):
 
    sol = open(sol_file,'w')
    sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(len(solution)) + '\n3 1 1 1\n')
    for i_point in range(len(solution)):
        sol.write(str(solution[i_point][0]) + " " + str(solution[i_point][1]) + " " + str(solution[i_point][2]) + "\n")
    sol.write('\nEnd\n')
    sol.close()

#: def write_sol()

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

#: def write_line()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-c", "--case",    dest="case",    default="None",
                      help="CASE", metavar="CASE")
                      
    (options, args)=parser.parse_args()

    main(options.case)
