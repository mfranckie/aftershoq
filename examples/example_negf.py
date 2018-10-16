'''
Created on 15 Feb 2018

@author: martin
'''

import sys
import os
import numpy as np

from aftershoq.structure import Structure
import aftershoq.structure.matpar as mp
from aftershoq.structure import Sgenerator, doping_modes as sgdm
from aftershoq.materials import *
from aftershoq.numerics.runplatf import Local
import aftershoq.utils.systemutil as su
import aftershoq.utils.debug as dbg
from aftershoq.interface import Inegf
from aftershoq.numerics import Paraopt, Gaussopt
from matplotlib import pyplot as pl
from aftershoq.qcls import EV2416


if __name__ == '__main__':
    
    # the working directory:
    path = "../../demo/EV2416/"
    path = os.getcwd()+"/"+path
    su.mkdir(path)
    
    # results directory:
    pathresults = "../../demo/results/"
    
    su.mkdir(pathresults)
    
    # Setup debugger:
    dbg.open(dbg.verb_modes['chatty'],pathresults + "/debug.log")
    dbg.debug("Debug file \n\n")
    dbg.flush()
    
    sep = '\n----------------------------------------------\n'
    
    print('--- Welcome to the "NEGF" demonstration of ---\n')
    print('               "AFTERSHOQ" \n\n')
    print('       Written by Martin Franckie 2018.')
    print('       Please give credit where credit')
    print('                   is due\n')
    print(sep)
    print('Creating semiconductor materials and alloys:\n')
    print('       ( All CBO relative to GaAs )\n')
    
    # create materials GaAs, AlAs, InAs:
    
    gaas = GaAs()
    
    for val in mp.valdict:
        print(val + " = " + str(gaas.params[mp.valdict[val]]))
        
    alas = AlAs()
    print("ALAs ELO = " + str( alas.params[mp.ELO] ))
    inas = InAs()
    print("\n")
    
    # create InGaAs/InAlAs lattice matched to InP:
    
    InGaAsLM = InGaAs()
    print(str(InGaAsLM) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(InGaAsLM.params[mp.valdict[val]]))
    
    InAlAsLM = AlInAs()
    print("\n" + str(InAlAsLM) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(InAlAsLM.params[mp.valdict[val]]))
    
    print("\nCBO = " + str(InAlAsLM.params[mp.Ec] - InGaAsLM.params[mp.Ec]) + "\n")
    
    # create Al_0.15Ga_0.85As 
    algaas = AlGaAs(x = 0.15)
    print("\n" + str(algaas) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(algaas.params[mp.valdict[val]]))
    
    print(sep)
    print('Creating a structure from generated materials:\n')
    print('[width, material, eta, lambda]')
    
    # creating a quantum-well structure
    '''
    s = Structure()
    
    # set interface roughness parameters for all interfaces:
    eta = 0.1
    lam = 2.0
    s.setIFR(eta, lam)
    
    # Add layers:
    s.addLayerWM(8.4,gaas)
    s.addLayerWM(3.1,algaas)
    s.addLayerWM(18.0,gaas)  # <--- this is layer 2!
    s.addLayerWM(1.8,algaas)
    s.addLayerWM(8.4,gaas)
    s.addLayerWM(3.1,algaas)
    s.addLayerWM(12.0,gaas)
    s.addLayerWM(1.8,algaas)
    
    # define doping layer
    zstart = 2; zend = 2.2; dopdens = 2e17; layer = 2
    
    # add a doping layer
    s.addDoping(zstart, zend, dopdens, layer)
    '''
    
    s = EV2416()
    s.addDoping(0, 1, 1e17, 1)
    
    print("Structure created : ")
    print(s)
    for l in s.layers:
        print(l.material.params[mp.Ec])
    
    print(sep)
    print('Generating N random structures:\n')
    N = eval(input('N = ?\n'))
    # define variations in composition, layer widths,
    # and doping location/density:
    dx = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    dw = [3,0,3,0,1,0]
    ddop = [0,0,1e16]
    
    # create a structure generator instance, based on our structure s above:
    sg = Sgenerator(s,dw,dx,ddop)
    
    sg.doping_mode = sgdm.get("FOLLOW LAYER VOLUME")
    
    # generate N random structures with the distribution in parameters defined above:
    #sg.genRanStructs(N)
    
    # generate N random structures, along the Hilbert curve with p = 5
    coords = sg.genRanHilbertStructs(N, 7)
    
    print ("Doping mode = " + str( sg.doping_mode) )
    
    for st in sg.structures:
        print(str(st.sid) + ' ' + str(st))
        print(st.dopings)
        print(st.getSheetDop())
        print("Doped layer:")
        li = st.layerIndex(st.dopings[0][0])
        z0 = st.layerPos(li)
        zf = z0 + st.layers[li].width
        print("Index=" + str(li) + ", (zi, zf) = (" + 
              str(z0) + ", " + str(zf) + ")")
        print()
        print("\n")
        
        
    print(sep)
    
    proceed = eval(input("Proceed to create directory tree?\nYes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
        dbg.close()
        exit()
    
    # set numerical parameters:
    
    # the binaries (change to your bin folder):
    binpath = "/Users/martin/git/NEGFT8/bin/"
    
    print('Creating directory tree in: ' + path + '.\n')
    print('Please change variable "binpath" to match your binaries!\n')
    
    # make the directory:
    su.mkdir(path)
    
    # Define the platform, local or cluster (e. g. Euler cluster at ETH Zuerich)
    pltfm = Local()
    #pltfm = Euler(1,"1:00")
    
    # negft interface:
    
    model = Inegf(binpath,pltfm,gaas)
    
    model.numpar["Nstates"] = 4
    model.numpar["NE"] = 200
    model.numpar["Nk"] = 200
    model.numpar["Tlattice"] = 300
    model.numpar["maxits"] = 30
    model.numpar["Nhist"] = 30
    model.numpar["gen"] = 5e-2
    model.numpar["Nh"] = 0
    model.numpar["omega0"] = 0.001
    model.numpar["domega"] = 0.030
    model.numpar["Nomega"] = 1
    model.numpar["efield0"] = 0.030
    model.numpar["defield"] = 0.001
    model.numpar["Nefield"] = 3
    model.numpar["Nper"] = 1
    model.numpar["fgr_omega0"] = 0
    model.numpar["fgr_omegaf"] = 0.030
    model.numpar["fgr_Nomega"] = 500
    
    # define the merit function as max gain/current density, from a dictionary of merit funcs.:
    #model.merit = model.merits.get("max gain")
    model.merit = model.merits.get("estimated gain")
    model.numpar['fgr_omega0'] = 0.010
    model.numpar['fgr_omegaf'] = 0.020
    
    print(sep + 'Starting simulations in directory tree? This will overwrite any previous results.')
    proceed = eval(input("Yes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
        dbg.close()
        exit()
    
    # execute the model for simulating the structures:
    model.datpath = '/Run/'
    model.runStructures(sg.structures, path)
    
    print("NEGF is running")
    model.waitforproc(1)
    
    print(sep + 'Gathering results to: ' + pathresults + '/results.log')
    
    # gather the results from the simulation of all structures:
    model.gatherResults(sg.structures, path, 
                        pathresults=pathresults, runprog=True)
    
    print(sep + 'First results acheived. Proceed with optimization?')
    proceed = eval(input("Yes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
        dbg.close()
        exit()
    
    # create optimization object
    tol, r, itmax, procmax = 0.001, 1.1, 20, 2
    #opt = Paraopt(tol,r,itmax,procmax)
    
    util = lambda mean, cov: -mean + np.sqrt(np.abs(np.diag(cov)))
    
    opt = Gaussopt(tol*sg.hutil.imax,itmax,procmax, sigma = 100, 
                   l = sg.hutil.imax*0.01, 
                   sigma_noise = 0.1, sigma_noise_max = 5,
                   padding = 20, utility=util)
    
    x0,y0 = opt.addEvaldPoints(model,sg,path,coords)
    
    print('Array of trial results:\n' + str(y0))
    print('evaluated at:\n' + str(x0))
    
    conv = opt.minimize(model, sg, path, pathresults)
    
    print("Optimization distances: \n" + str(opt.x))
    print("Optimization values: \n"    + str(opt.y))
    
    
    pl.plot(opt.x,opt.y)
    pl.show()
    
    print(sep + 'Program completed! Good bye!')
    dbg.close()
    
