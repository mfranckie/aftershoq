'''
Created on 15 Feb 2018

@author: martin
'''

from classes import Structure, Sgenerator, MaterialPar as mp
from interface import Inegf, NumPar as np
from qclutil import MaterialUtil as mu
from runplatf import Local
from systemutil import SystemUtil as su
import os

if __name__ == '__main__':
    
    
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
    
    GaAs = mu.createGaAs()
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
    x = 0.15
    AlGaAs = mu.createAlGaAs(x, "Al0.15Ga0.85As")
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
    s.addLayerMW(8.4,GaAs)
    s.addLayerMW(3.1,AlGaAs)
    s.addLayerMW(18.0,GaAs)  # <--- this is layer 2!
    s.addLayerMW(1.8,AlGaAs)
    
    # define doping layer
    zstart = 2; zend = 2.2; dopdens = 2e17; layer = 2
    
    # add a doping layer
    s.addDoping(zstart, zend, dopdens, layer)
    
    print s
    
    print sep
    print 'Generating N random structures:\n'
    N = input('N = ?')
    
    # define variations in composition (not implemented yet), layer widths,
    # and doping location/density:
    dx = [0.05, 0.05, 0.05, 0.05]
    dw = [10,10,10,10]
    ddop = [0.1,1,1e17]
    
    # create a structure generator instance, based on our structure s above:
    sg = Sgenerator(s,dw,dx,ddop)
    
    # generate N random structures with the distribution in parameters defined above:
    sg.genRanStructs(N)
    for st in sg.structures:
        print str(st.sid) + ' ' + str(st)
        
    print sep
    
    proceed = input("Proceed to create directory tree?\nYes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        exit()
    
    # set numerical parameters:
    
    # the binaries (change to your bin folder):
    binpath = "/Users/martin/git/NEGFT8/bin/"
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
    numpar[np.Nnu] = 3
    numpar[np.Nper] = 2
    
    # Define the platform, local or cluster (e. g. Euler cluster at ETH Zuerich)
    pltfm = Local()
    #pltfm = Euler(6,"1:00")
    
    # select model for transport and em raditaion simulations:
    
    # NEGF interface:
    model = Inegf(binpath,pltfm,numpar,GaAs)
    
    # define the merit function as max gain/current density, from a dictionary of merit funcs.:
    model.merit = Inegf.merits.get("(max gain)/(current density)")
    
    print sep + 'Starting simulations in directory tree?'
    proceed = input("Yes = 1\nExit = 0")
    if proceed == 0:
        print sep + "User exit. Good bye!"
        exit()
    # execute the model for simulating the structures:
    model.runStructures(sg.structures, path)
    
    model.waitforproc(5, "negft running....")
    
    print sep + 'Gathering results to: ' + path + 'results.log'
    
    # gather the results from the simulation of all structures:
    model.gatherResults(sg.structures, path)
    
    print sep + 'Program completed! Good bye!'
    
    