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
from structure.materials import GaAs, AlGaAs
from numerics.runplatf import Local
from utils.systemutil import SystemUtil as su
from utils.debug import Debugger as dbg
from interface.isewlab import Isewlab
from numerics.paraopt import Paraopt
from matplotlib import pyplot as pl
import numpy as np


if __name__ == '__main__':
    
    # the binaries (change to your bin folder):
    binpath = "/usr/local/bin/sewlab_MAC"
    
    # the working directory:
    path = "../../demo/"
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
    
    # create materials:
    # create GaAs:
    gaas = GaAs()
    
    
    # create Al_0.15Ga_0.85As 
    algaas = AlGaAs(x = 0.15)
    print "\n" + str(algaas) + ":\n"
    for val in mp.valdict:
        print val + " = " + str(algaas.params[mp.valdict[val]])
    
    print sep
    print 'Creating a structure from generated materials:\n'
    print '[width, material, eta, lambda]'
    
    # creating a two quantum-well structure
    mystruct = Structure()
    
    # set interface roughness parameters for all interfaces:
    eta = 0.1
    lam = 2.0
    mystruct.setIFR(eta, lam)
    
    # Add layers:
    mystruct.addLayerMW(5.1,algaas)
    mystruct.addLayerMW(9.3,gaas)
    mystruct.addLayerMW(1.0,algaas)
    mystruct.addLayerMW(10.5,gaas)
    mystruct.addLayerMW(3.3,algaas)
    mystruct.addLayerMW(8.7,gaas)
    mystruct.addLayerMW(4.2,algaas)
    mystruct.addLayerMW(16.5,gaas)

    # when using sewlab, doping will cover entire layer.
    # to make a different doping thickness, create a 
    # separate doping layer.
    mystruct.addDoping(0, 16.5, 2e16, 7)
    
    print "Structure created : "
    print mystruct
    
    print sep
    
    print 'Creating directory tree in: ' + path + '.\n'
    print 'Please change variable "binpath" to match your binaries!\n'
    
    # Define the platform, local or cluster (e. g. Euler cluster at ETH Zuerich)
    pltfm = Local()
    #pltfm = Euler(1,"1:00")
    
    # sewlab interface: (add parameter version '4.6.X' if you are not running the default 4.6.4
    model = Isewlab(binpath,pltfm, gaas)
    
    # to change parameters, change the dictionaries in Isewlab:
    
    # determines the ammount of output produced (also: "Verbose")
    model.numpar["verbosity"] = "Silent"
    
    # The electric field in the simulation
    efield0 = -1
    defield = -2
    Nefield = 5
    model.numpar["efield0"] = efield0
    
    # Set the elecron and lattice temperatures
    Te = 100 # Kelvin
    Tl = 77  # Kelvin
    model.setTe(Te)         # electron temperature (set to 1500 K if kinetic balance is used)
    model.setTlattice(Tl)   # lattice tempreature (~cryostat temperature)
    
    model.useKinBal(False) # Use of kinetic balance?
    model.useSuperself(False) # Use of superself?
    model.computeLight(False) # Use of light computations (power, photon-driven current)
    
    # to create input files with default parameters (automatically done via runStructures() below):
    model.writeSampleFile(mystruct, path)
    model.writeScriptFile(path)
    
    print sep + 'Starting simulations in directory tree (sewlab version ' + model.version + ' ). This will overwrite any previous results.'
    
    # for looping over bias points:
    for i in range(0,Nefield):
        model.numpar["efield0"] = efield0 + defield*i
    
        # execute the model for simulating the structures and wait:
        # (Comment out these next two lines if you just want to gather results which were already simulated)
        model.runStructures( [mystruct], path )
        model.waitforproc(5, "sewlab is running.... iteration = " + str( i ))
    
    print sep + 'Gathering results to: ' + path + '/results.log'
    
    # gather the results from the simulation of all structures:
    model.gatherResults( [mystruct], path )
    
    print sep + 'Results acheived!'
    
    # now, the results for each bias are stored in mystruct.results, .dipoles, .energies, .populations, and .rates
    results = np.array( mystruct.results )
    descr = ['E (kV/cm)', 'j (A/cm^2)', 'Max gain (1/cm)', 'Peak energy. (meV)']
    print descr
    print results
    # plots the I-V:
    pl.plot(-results[:,0], results[:,1],'-*')
    pl.xlabel("Electric field (kV/cm)")
    pl.ylabel("Current density (A/cm^2)")
    pl.show()
    
    
    print sep + 'Program completed! Good bye!'
    dbg.close()
    