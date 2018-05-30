'''
Created on 15 Feb 2018

@author: martin
'''

import os
import sys
# change path as apropriate
path_to_aftershoq = os.getcwd()
sys.path.append(path_to_aftershoq)
sys.path.append(path_to_aftershoq + '/hilbert_curve/')

from structure.classes import Structure, MaterialPar as mp
from structure.materials import *
from structure.sgenerator import Sgenerator
from interface.inegf import Inegf, NumPar as np
from utils.qclutil import MaterialUtil as mu
from numerics.runplatf import Local
from utils.systemutil import SystemUtil as su
from utils.debug import Debugger as dbg
from interface.isewself import Isewself
from numerics.paraopt import Paraopt
import numpy




if __name__ == '__main__':
    
    # Setup debugger:
    
    dbg.open(dbg.verb_modes['verbose'],"debug.log")
    dbg.debug("Debug file\n\n")
    
    sep = '\n----------------------------------------------\n'
    
    print '------ Welcome to the demonstration of -------\n'
    print '               "AFTERSHOQ" \n\n'
    print '       Written by Martin Franckie 2018.'
    print '       Please give credit where credit'
    print '                   is due\n'
    print sep
    print 'Creating semiconductor materials and alloys:\n'
    print '       ( All CBO relative to GaAs )\n'
    
    # create materials GaAs, AlAs, InAs:
    
    gaas = GaAs()
    print str(gaas) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(gaas.params[mp.valdict[val]])
    
    alas = AlAs()
    print "\n" + str(alas) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(alas.params[mp.valdict[val]])
        
    inas = InAs()
    print "\n" + str(inas) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(inas.params[mp.valdict[val]])
    
    # create InGaAs/InAlAs lattice matched to InP:
    
    x = 0.47
    lmname = "In" + str(x) + "Ga" + str(1-x) + "As"
    InGaAsLM = InGaAs(lmname, x)
    print "\n" + str(InGaAsLM) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(InGaAsLM.params[mp.valdict[val]])
    
    x = 0.48
    lmname = "In" + str(x) + "Al" + str(1-x) + "As"
    InAlAsLM = InAlAs(lmname,x)
    print "\n" + str(InAlAsLM) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(InAlAsLM.params[mp.valdict[val]])
    
    # create Al_0.15Ga_0.85As 
    algaas = AlGaAs(x = 0.15)
    print "\n" + str(algaas) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(algaas.params[mp.valdict[val]])
    
    print sep
    print 'Creating a structure from generated materials:\n'
    print '[width, material, eta, lambda]'
    
    # creating a two quantum-well structure
    s = Structure()
    
    # set interface roughness parameters for all interfaces:
    eta = 0.1
    lam = 2.0
    s.setIFR(eta, lam)
    
    # Add layers:
    s.addLayerMW(8.4,gaas)
    s.addLayerMW(3.1,algaas)
    s.addLayerMW(18.0,gaas)  # <--- this is layer 2!
    s.addLayerMW(1.8,algaas)
    s.addLayerMW(8.4,gaas)
    s.addLayerMW(3.1,algaas)
    s.addLayerMW(18.0,gaas)
    s.addLayerMW(1.8,algaas)
    
    # define doping layer
    zstart = 2; zend = 2.2; dopdens = 2e17; layer = 2
    
    # add a doping layer
    s.addDoping(zstart, zend, dopdens, layer)
    
    print s
    
    print sep
    print 'Generating N random structures:\n'
    N = input('N = ?\n')
    
    # define variations in composition (not implemented yet), layer widths,
    # and doping location/density:
    dx = [0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1]
    dw = [1,1,1,1,1,1,1,1]
    ddop = [0,0,0]
    
    # create a structure generator instance, based on our structure s above:
    sg = Sgenerator(s,dw,dx,ddop)
    
    # generate N random structures with the distribution in parameters defined above:
    #sg.genRanStructs(N)
    
    # generate N random structures, along the Hilbert curve with p = 5
    coords = sg.genRanHilbertStructs(N, 5)
    
    for st in sg.structures:
        print str(st.sid) + ' ' + str(st)
        
    print sep
    
    proceed = input("Proceed to create directory tree?\nYes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        dbg.close()
        exit()
    
    # set numerical parameters:
    
    # the binaries (change to your bin folder):
    # NEGF:
    #binpath = "/Users/martin/git/NEGFT8/bin/"
    
    # SEWSELF
    binpath = "/Users/martin/git/Sewself_JF/bin/"
    
    # the execution directory:
    path = "../demo"
    path = os.getcwd()+"/"+path
    
    print 'Creating directory tree in: ' + path + '.\n'
    print 'Please change variable "binpath" to match your binaries!\n'
    
    # make the directory:
    su.mkdir(path)
    
    # list of numerical parameters; to be expanded for future model interfaces
    numpar = []
    np.initList(numpar)
    np.setDefault(numpar)
    
    # to change numerical paramters, use the following syntax:
    numpar[np.Nnu] = 1
    numpar[np.Nper] = 1
    numpar[np.Niter] = 1
    numpar[np.Nhist] = 1
    
    # Define the platform, local or cluster (e. g. Euler cluster at ETH Zuerich)
    pltfm = Local()
    #pltfm = Euler(6,"1:00")
    
    # select model for transport and em raditaion simulations:
    
    # NEGF interface:
    #model = Inegf(binpath,pltfm,numpar,GaAs)
    
    # sewself interface:
    material_list = [GaAs,AlGaAs]
    model = Isewself(binpath,pltfm,material_list)
    model.writeStructFile(s, path)
    model.writeParameterFile(path)
    model.writeSewselfPar(path)
    
    # define the merit function as max gain/current density, from a dictionary of merit funcs.:
    model.merit = model.merits.get("(max gain)/(current density)")
    
    print sep + 'Starting simulations in directory tree?'
    proceed = input("Yes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        dbg.close()
        exit()
        
    # execute the model for simulating the structures:
    model.runStructures(sg.structures, path)
    
    model.waitforproc(1, str(model) + " running....")
    
    print sep + 'Gathering results to: ' + path + '/results.log'
    
    # gather the results from the simulation of all structures:
    model.gatherResults(sg.structures, path)
    
    print sep + 'First results acheived. Proceed with optimization of gain?'
    proceed = input("Yes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        dbg.close()
        exit()
    
    # collect results from trial points
    x0 = []
    [x0.append( sg.hutil.interp_dist_from_coords( c ) ) for c in coords]
    x0 = numpy.array(x0)
    x0.sort()
    x0 = x0.tolist()
    print x0
    y0 = []
    xi = 0
    for i in range(0,len(x0)):
        try:
            y0.append( -float(model.getMerit(sg.structures[i],path)) )
        except( ValueError ):
            del x0[xi]
            xi-=1
        xi +=1
    
    print 'Array of trial results:\n' + str(y0)

    # create optimization object
    tol, r, itmax, procmax = 1, 1.3, 10, 2 
    opt = Paraopt(tol,r,itmax,procmax,x0,y0)
    
    conv = opt.minimize(model, sg, path)
    
    
    
    print sep + 'Program completed! Good bye!'
    dbg.close()
    
    