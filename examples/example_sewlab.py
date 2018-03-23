'''
Created on 15 Feb 2018

@author: martin
'''
import sys
import os
import numpy
# change path as apropriate
path_to_aftershoq = os.getcwd()
sys.path.append(path_to_aftershoq)
sys.path.append(path_to_aftershoq + '/hilbert_curve/')

from structure.classes import Structure, MaterialPar as mp
from structure.sgenerator import Sgenerator
from utils.qclutil import MaterialUtil as mu
from numerics.runplatf import Local
from utils.systemutil import SystemUtil as su
from utils.debug import Debugger as dbg
from interface.isewlab import Isewlab
from numerics.paraopt import Paraopt
from matplotlib import pyplot as pl


if __name__ == '__main__':
    
    # the working directory:
    path = "../../demo"
    path = os.getcwd()+"/"+path
    su.mkdir(path)
    
    # Setup debugger:
    dbg.open(dbg.verb_modes['verbose'],path+"/debug.log")
    dbg.debug("Debug file \n\n")
    dbg.flush()
    
    sep = '\n----------------------------------------------\n'
    
    print '--- Welcome to the "sewlab" demonstration of ---\n'
    print '               "AFTERSHOQ" \n\n'
    print '       Written by Martin Franckie 2018.'
    print '       Please give credit where credit'
    print '                   is due\n'
    print sep
    print 'Creating semiconductor materials and alloys:\n'
    print '       ( All CBO relative to GaAs )\n'
    
    # create materials GaAs, AlAs, InAs:
    
    GaAs = mu.createGaAs()
    GaAs.params[mp.Eg] = 1.16
    AlAs = mu.createAlAs()
    InAs = mu.createInAs()
    
    # create InGaAs/InAlAs lattice matched to InP:
    
    x = 0.47
    lmname = "In" + str(x) + "Ga" + str(1-x) + "As"
    InGaAsLM = mu.createGaInAs(x,lmname)
    print str(InGaAsLM) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(InGaAsLM.params[mp.valdict[val]])
    
    x = 0.48
    lmname = "In" + str(x) + "Al" + str(1-x) + "As"
    InAlAs = mu.createAlInAs(x,lmname)
    print "\n" + str(InAlAs) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(InAlAs.params[mp.valdict[val]])
    
    # create Al_0.15Ga_0.85As 
    x = 0.20
    AlGaAs = mu.createAlGaAs(x)
    AlGaAs.params[mp.Eg] = 1.90
    print "\n" + str(AlGaAs) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(AlGaAs.params[mp.valdict[val]])
    
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
    s.addLayerMW(5.1,AlGaAs)
    s.addLayerMW(9.3,GaAs)
    s.addLayerMW(1.0,AlGaAs)
    s.addLayerMW(10.5,GaAs)
    s.addLayerMW(3.3,AlGaAs)
    s.addLayerMW(8.7,GaAs)
    s.addLayerMW(4.2,AlGaAs)
    s.addLayerMW(16.5,GaAs)

    s.addDoping(0, 16.5, 2e16, 7)
    
    print "Structure created : "
    print s
    
    print sep
    print 'Generating N random structures:\n'
    N = input('N = ?\n')
    
    # define variations in composition, layer widths,
    # and doping location/density:
    dx = [0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0]
    dw = [0,0,2,0,2,0,0,0]
    ddop = [0,0,0]
    
    # create a structure generator instance, based on our structure s above:
    sg = Sgenerator(s,dw,dx,ddop)
    
    # generate N random structures with the distribution in parameters defined above:
    #sg.genRanStructs(N)
    
    # generate N random structures, along the Hilbert curve with p = 5
    coords = sg.genRanHilbertStructs(N, 7)
    
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
    binpath = "/usr/local/bin/"
    
    print 'Creating directory tree in: ' + path + '.\n'
    print 'Please change variable "binpath" to match your binaries!\n'
    
    # make the directory:
    su.mkdir(path)
    
    # Define the platform, local or cluster (e. g. Euler cluster at ETH Zuerich)
    pltfm = Local()
    #pltfm = Euler(1,"1:00")
    
    # sewself interface:
    material_list = [GaAs,AlGaAs]
    model = Isewlab(binpath,pltfm)
    model.numpar["verbosity"] = "Silent"
    model.numpar["efield0"] = -11
    
    # to change parameters, change the dictionaries in isewlab:
        
    
    # to create input files with default parameters:
    
    
    # define the merit function as max gain/current density, from a dictionary of merit funcs.:
    
    
    print sep + 'Starting simulations in directory tree? This will overwrite any previous results.'
    proceed = input("Yes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        dbg.close()
        exit()
    
    # execute the model for simulating the structures:
    model.runStructures(sg.structures, path)
    
    model.waitforproc(5, "sewlab is running....")
    
    print sep + 'Gathering results to: ' + path + '/results.log'
    
    # gather the results from the simulation of all structures:
    model.gatherResults(sg.structures, path)
    
    model.merit = model.merits.get("max gain")
    print model.merit
    print model.merits
    model.target = 0.20
    
    print sep + 'First results acheived. Proceed with optimization?'
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
    tol, r, itmax, procmax = 18, 1.1, 100, 2
    
    opt = Paraopt(tol,r,itmax,procmax,x0,y0)
    
    conv = opt.minimize(model, sg, path)
    
    print "Optimization distances: \n" + str(opt.x)
    print "Optimization values: \n"    + str(opt.y)
    
    pl.plot(opt.x,opt.y)
    pl.show()
    
    print sep + 'Program completed! Good bye!'
    dbg.close()
    