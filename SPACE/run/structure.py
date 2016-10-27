#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import sys, os, copy

sys.path.append(os.environ['STRUCTURE_RUN'])

import numpy, argparse
from mpi4py import MPI
from baseclasses import *
from tacs import *
from pyOpt import *

from .. import io   as spaceio
from .. import util as spaceutil
from interface import INT       as SPACE_INT
from load      import load      as spaceload

# ----------------------------------------------------------------------
#  Structure Simulation
# ----------------------------------------------------------------------

def structure(config):

    gcomm = comm = MPI.COMM_WORLD

    # local copy
    konfig = copy.deepcopy(config)

    # Load
    spaceutil.surf2pressuresol(konfig)
    SPACE_INT(konfig)
    info = spaceload(konfig)

    # Material properties

    rho_aluminum = 2780
    E_aluminum = 72e9
    ys_aluminum = 324e6
    nu_aluminum = 0.33
    kcorr_aluminum = 5.0/6.0

    rho_titanium = 4460
    E_titanium = 110e9
    ys_titanium = 869e6
    nu_titanium = 0.31
    kcorr_titanium = 5.0/6.0

    t = 0.01
    tMin = 0.0016 # 0.0016 # 1/16"
    tMax = 0.05 # 0.020

    nx = float(konfig.ACCELERATION_X)
    ny = float(konfig.ACCELERATION_Y)
    nz = float(konfig.ACCELERATION_Z)

    loadFactor = numpy.sqrt(nx*nx+ny*ny+nz*nz)
    gravityVector = -9.81 * numpy.array([-nx,-nz,-ny])/numpy.sqrt(nx*nx+ny*ny+nz*nz) # Change of Frame: to Structure Frame



    KSWeight = 80.0
    SPs = [StructProblem('lc0', loadFactor=loadFactor, loadFile=konfig.LOAD_FILENAME, evalFuncs=['mass','ks0'])]
    #SPs = [StructProblem('lc0', loadFactor=loadFactor, loadFile=konfig.LOAD_FILENAME, evalFuncs=['mass','ks0','ks1','ks2'])]
    numLoadCases = len(SPs)

    structOptions = {'transferSize':0.5, 'transferGaussOrder':3}

    FEASolver = pytacs.pyTACS(konfig.STRUCT_SURFACE + '.bdf', comm=comm, options=structOptions)

    SKIN_FUSE_U = ['FUSE:TOP','FUSE:LFT','CTAIL:LOW','CTAIL_T','FUSE_R']
    SKIN_FUSE_L = ['FUSE:BOT','FLAP:UPP','FLAP:LOW','FLAP_T','FUSE_F']
    SKIN_WING_U = ['LWING:UPP','LWING_T::0']
    SKIN_WING_L = ['LWING:LOW','LWING_T::1']
    SKINS = SKIN_FUSE_U + SKIN_FUSE_L + SKIN_WING_U + SKIN_WING_L + ['MSKINC:a','MSKINC:b']

    JUNCTIONS = ['FLAP_FUSE','LWING_FUSE','CTAIL_FUSE']

    FRAMES = ['MFRAME:00','MFRAME:01','MFRAME:02','MFRAME:03','MFRAME:04','MFRAME:05','MFRAME:06','MFRAME:07','MFRAME:08','MFRAME:09',
        'MFRAME:10','MFRAME:11','MFRAME:12']
    LONGERONS = ['MLONG:02:2','MLONG:00:3','MLONG:01:3','MLONG:02:3','MLONG:00:4','MLONG:01:4']
    RIBS = ['MRIBF:00','MRIBF:01','MRIBF:02','MRIBF:03','MRIBF:04','MRIBF:05','MRIBF:06','MRIBF:07',
        'MRIBV:00','MRIBV:01','MRIBV:02','MRIBV:03','MRIBV:04','MRIBV:05','MRIBV:06','MRIBV:07','MRIBV:08','MRIBV:09',
        'MRIBW:00','MRIBW:01','MRIBW:02','MRIBW:03','MRIBW:04','MRIBW:05']
    SPARS = ['MSPARF:02','MSPARF:03',
        'MSPARV:00','MSPARV:01',
        'MSPARC:00','MSPARC:06','MSPARC:07',
        'MSPARW:00','MSPARW:02','MSPARW:08','MSPARW:09']
    STRINGERS = ['MSTRINGC:01','MSTRINGC:02','MSTRINGC:03','MSTRINGC:04','MSTRINGC:05',
        'MSTRINGW:01','MSTRINGW:03','MSTRINGW:04','MSTRINGW:05','MSTRINGW:06','MSTRINGW:07']
    MEMBERS = FRAMES + LONGERONS + RIBS + SPARS + STRINGERS

    corresp = [-1 for index in range(len(FEASolver.selectCompIDs(include=SKINS+JUNCTIONS+MEMBERS)[0]))]
    dv = -1;
    SKIN_IDS = FEASolver.selectCompIDs(include=SKINS)[0]
    print SKIN_IDS
    for i in range(len(SKIN_IDS)):
        dv_name = "SKIN_" + str(i)
        FEASolver.addDVGroup(dv_name, include = SKIN_IDS[i])
        dv = dv+1;
        corresp[SKIN_IDS[i]] = dv;
    JUNCTION_IDS = FEASolver.selectCompIDs(include=JUNCTIONS)[0]
    for i in range(len(JUNCTION_IDS)):
        dv_name = "JUNCTION_" + str(i)
        FEASolver.addDVGroup(dv_name, include = JUNCTION_IDS[i])
        dv = dv+1;
        corresp[JUNCTION_IDS[i]] = dv;
    for i in range(len(MEMBERS)):
        dv_name = MEMBERS[i]
        FEASolver.addDVGroup(dv_name, include = MEMBERS[i])
        dv = dv+1;
        MEMBERS_IDS_I = FEASolver.selectCompIDs(include=MEMBERS[i])[0]
        for k in range(len(MEMBERS_IDS_I)):
            corresp[MEMBERS_IDS_I[k]] = dv;

    #print(corresp)
    #print(len(FEASolver.selectCompIDs(include=SKINS+JUNCTIONS+MEMBERS)[0]))

    def conCallBack(dvNum, compDescripts, userDescript, specialDVs, **kargs):
        con = constitutive.isoFSDTStiffness(rho_aluminum, E_aluminum, nu_aluminum, kcorr_aluminum, ys_aluminum, t, dvNum, tMin, tMax)
        # if userDescript in ['JUNCTIONS','FRAMES','LONGERONS','SPARS']:
        #     con = constitutive.isoFSDTStiffness(rho_2024, E_2024, nu, kcorr, ys_2024, 0.05, dvNum, 0.0016, 0.05)
        scale = [100.0]
        return con, scale

    FEASolver.createTACSAssembler(conCallBack)

    # Mass Functions
    FEASolver.addFunction('mass', functions.StructuralMass)

    # KS Functions
    ks0 = FEASolver.addFunction('ks0', functions.AverageKSFailure, KSWeight=KSWeight, loadFactor=5.0)
    #ks0 = FEASolver.addFunction('ks0', functions.AverageKSFailure, KSWeight=KSWeight, include=RIBS+SPARS+FRAMES+LONGERONS+WING_BOX, loadFactor=loadFactor)
    #ks1 = FEASolver.addFunction('ks1', functions.AverageKSFailure, KSWeight=KSWeight, include=SKIN_U+STRINGERS_U, loadFactor=loadFactor)
    #ks2 = FEASolver.addFunction('ks2', functions.AverageKSFailure, KSWeight=KSWeight, include=SKIN_L+STRINGERS_L, loadFactor=loadFactor)

    # Load Factor
    FEASolver.setOption('gravityVector',gravityVector.tolist())
    for i in range(numLoadCases):
       FEASolver.addInertialLoad(SPs[i])

    # Optimize

    def obj(x):
        '''Evaluate the objective and constraints'''
        funcs = {}
        FEASolver.setDesignVars(x)
        for i in range(numLoadCases):
            FEASolver(SPs[i])
            FEASolver.evalFunctions(SPs[i], funcs)
        if comm.rank == 0:
            print funcs
        f = funcs['lc0_mass']
        g = []
        g.append(funcs['lc0_ks0'])
        #g.append(funcs['lc0_ks1'])
        #g.append(funcs['lc0_ks2'])
        fail = 0
        return f, g, fail

    def sens(x, f, g):
        '''Evaluate the objective and constraint sensitivities'''
        funcsSens = {}
        
        x_cur = numpy.zeros(len(x))
        FEASolver.structure.getDesignVars(x_cur)
        if not numpy.array_equal(x_cur,x):
            print 'need recompute'
            FEASolver.setDesignVars(x)
            for i in range(numLoadCases):
                FEASolver(SPs[i])

        for i in range(numLoadCases):
            FEASolver.evalFunctionsSens(SPs[i], funcsSens)

        df1 = funcsSens['lc0_mass'][FEASolver.varSet]

        dg1 = []
        dg1.append(funcsSens['lc0_ks0'][FEASolver.varSet])
        #dg1.append(funcsSens['lc0_ks1'][FEASolver.varSet])
        #dg1.append(funcsSens['lc0_ks2'][FEASolver.varSet])

        df = [0.0]*len(df1)
        for i in range(0,len(df1)):
            df[i] = df1[i] # NEEDED ?????????????????

        dg = numpy.zeros([len(dg1),len(dg1[0])])
        for i in range(0,len(dg1[0])):
            dg[0][i] = dg1[0][i]
            #dg[1][i] = dg1[1][i]
            #dg[2][i] = dg1[2][i]

        fail = 0

        return df, dg, fail


    optimize = True
    if optimize:

        # Set up the optimization problem
        optProb = Optimization('Mass min', obj)
        optProb.addObj('lc0_mass')
        FEASolver.addVariablesPyOpt(optProb)

        for i in range(numLoadCases):
            for j in xrange(1):   ###################################################
                optProb.addCon('%s_ks%d'% (SPs[i].name, j), upper=1.0)

        if comm.rank == 0:
            print optProb
        #optProb.printSparsity()

        # Instantiate Optimizer (SNOPT) & Solve Problem

        snopt = SNOPT()
        #snopt.setOption('Major iterations limit',3)
        [obj_fun, x_dvs, inform] = snopt(optProb, sens_type=sens)
        #print optProb.solution(0)
        optProb.write2file(outfile='pySNOPT.txt', disp_sols=False, solutions=[0])

        # Setting the final solution

        FEASolver.setDesignVars(x_dvs) # NEEDED ?????????????????

    # Getting the final solution

    for i in range(numLoadCases):
        FEASolver(SPs[i]) # NEEDED ?????????????????

    dvs_file = 'dvs.dat'
    dvs = open(dvs_file,'w')
    for i in range(len(corresp)):
        dvs.write('%f\n' % x_dvs[corresp[i]])
    dvs.close()

    # funcs = {}
    # FEASolver.evalFunctions(SPs[0], funcs)

    FEASolver.writeBDFForces(SPs[0], "visualize_forces.bdf")
    FEASolver.writeSolution()
    # #FEASolver.writeBDF("crm_wing_design_tacs_final.bdf")


    # info out
    info = spaceio.State(info)
    #info.FILES.MASS = konfig.MASS
    return info

#: def structure()

