'''
Created on 15 Feb 2018

@author: martin
'''
import sys
import os
import numpy

from aftershoq.structure import Structure
import aftershoq.structure.matpar as mp
from aftershoq.materials import *
from aftershoq.structure import Sgenerator
from aftershoq.numerics.runplatf import Local
import aftershoq.utils.systemutil as su
import aftershoq.utils.debug as dbg
from aftershoq.interface import Isewself, Inegf
from aftershoq.numerics import Paraopt, Gaussopt
from matplotlib import pyplot as pl
import aftershoq.utils.const as const


if __name__ == '__main__':
    
    # the working directory:
    path = "/Users/martin/Projects/Optimization/ASQW/ASQWchi2/Mehrans/v3/"
    #path = os.getcwd()+"/"+path
    su.mkdir(path)
    
    # Setup debugger:
    dbg.open(dbg.verb_modes['verbose'],path+"/debug.log")
    dbg.debug("Debug file \n\n")
    dbg.flush()
    
    sep = '\n----------------------------------------------\n'
    
    print('--- Welcome to the "sewself" demonstration of ---\n')
    print('               "AFTERSHOQ" \n\n')
    print('       Written by Martin Franckie 2018.')
    print('       Please give credit where credit')
    print('                   is due\n')
    print(sep)
    print('Creating semiconductor materials and alloys:\n')
    print('       ( All CBO relative to GaAs )\n')
    
    # create materials GaAs and two AlGaAs:
    
    gaas = GaAs()
    alas = AlAs()
    algaas1 = AlGaAs(name = "AlGaAs1", x = 0.22)
    algaas2 = AlGaAs(name = "AlGaAs2", x = 0.43)
    
    print(sep)
    print('Creating a structure from generated materials:\n')
    print('[width, material, eta, lambda]')
    
    # creating a two quantum-well structure
    s = Structure()
    
    # set interface roughness parameters for all interfaces:
    eta = 0.1
    lam = 2.0
    s.setIFR(eta, lam)
    
    # Add layers:
    s.addLayerWM(7.4, algaas2)
    s.addLayerWM(2.5, gaas) 
    s.addLayerWM(0.5, gaas) # <-- doped layer 2
    s.addLayerWM(2.5, gaas) 
    s.addLayerWM(8.7, algaas1)
    s.addLayerWM(7.4,algaas2)
    
    # define doping layer
    zstart = 0.01; zend = 0.499; dopdens = 2e19; layer = 2
    
    # add a doping layer
    s.addDoping(zstart, zend, dopdens, layer)
    
    print("Structure created : ")
    print(s)
    
    print(sep)
    print('Generating N random structures:\n')
    N = eval(input('N = ?\n'))
    
    # define variations in composition, layer widths,
    # and doping location/density:
    dx = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    dw = [0, 2, 0, 2, 4, 0]
    ddop = [0,0,0]
    
    # create a structure generator instance, based on our structure s above:
    sg = Sgenerator(s,dw,dx,ddop)
    
    # generate N random structures with the distribution in parameters defined above:
    #sg.genRanStructs(N)
    
    # generate N random structures, along the Hilbert curve with p = 7
    
    coords = sg.genRanHilbertStructs(N, 7)
    
    for st in sg.structures:
        print(str(st.sid) + ' ' + str(st))
        print(st.dopings)
        
    print(sep)
    
    proceed = eval(input("Proceed to create directory tree?\nYes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
        dbg.close()
        exit()
    
    # set numerical parameters:
    
    # the binaries (change to your bin folder):
    binpath = "/Users/martin/git/Sewself_JF/bin/"
    pathnegf = "/Users/martin/git/NEGFT8/bin/"
    
    print('Creating directory tree in: ' + path + '.\n')
    print('Please change variable "binpath" to match your binaries!\n')
    
    # make the directory:
    su.mkdir(path)
    
    # Define the platform, local or cluster (e. g. Euler cluster at ETH Zuerich)
    pltfm = Local()
    #pltfm = Euler(1,"1:00")
    
    # sewself interface:
    material_list = [gaas,alas]
    commands = ['e', 'd', 'a210', 'a310', 'a320']
    
    
    
    model = Isewself(binpath,pltfm,material_list,commands = commands)
    #model = Inegf(pathnegf,pltfm,gaas)
    
    #model.writeWannier(s, path)
    
    
    E1 = 0.162
    E2 = 0.1536
    gamma = 0.01
    
    model.target = [E1, E2, gamma]
    
    
    
    
    
    #model.numpar["efield0"] = 0.001
    #model.numpar["Nomega"] = 1
    #model.numpar["gen"] = 1e-3
    
    # to change parameters, change the dictionaries in isewself:
    model.sewselfpar["emin"] = -0.0001
    model.sewselfpar["emax"] = 0.34
    model.numpar["efield"] = -0.0
    model.sewselfpar["Espan"] = 1
    # (material parameters are atuo-generated from materials in material_list)
    
    # define the merit function as max gain/current density, from a dictionary of merit funcs.:
    model.merit = model.merits.get("Chi2")
    
    
    model.numpar['lattice_temp'] = 77
    model.numpar['el_temp'] = 77
    
    print(sep + 'Starting simulations in directory tree? This will overwrite any previous results.')
    proceed = eval(input("Yes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
        dbg.close()
        exit()
    
    #model.runStructures([s], path)
    #model.waitforproc(2, "Simulation is running for orig....")
    
    # execute the model for simulating the structures:
    model.runStructures([s] + sg.structures, path)
    
    model.waitforproc(2, "Simulation is running....")
    
    print(sep + 'Gathering results to: ' + path + '/results.log')
    
    # gather the results from the simulation of all structures:
    model.gatherResults([s] + sg.structures, path, runprog=True)
    
   
    
    
    
    print(sep + 'First results acheived. Proceed with optimization?')
    proceed = eval(input("Yes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
        dbg.close()
        exit()
        
    # create optimization object
    print("imax = " + str(sg.hutil.imax))
    
    
    tol, r, itmax, procmax = 0.0001*sg.hutil.imax, 1.1, 100, 5
    
    opt = Paraopt(tol,r,itmax,procmax)
    
    opt.addEvaldPoints(model, sg, path, coords)
    
    conv = opt.minimize(model, sg, path)
    
    
    print(sep + 'Program completed! Good bye!')
    dbg.close()
    
