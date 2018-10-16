'''
Created on 15 Feb 2018

@author: martin
'''

import os

from aftershoq.structure import Structure
import aftershoq.structure.matpar as mp
from aftershoq.materials import *
from aftershoq.qcls import *
from aftershoq.structure import Sgenerator
from aftershoq.interface import Inegf
from aftershoq.numerics.runplatf import Local
import aftershoq.utils.systemutil as su
import aftershoq.utils.debug as dbg
from aftershoq.interface import Isewself
from aftershoq.numerics import Paraopt
import numpy




if __name__ == '__main__':
    
    # Setup debugger:
    
    dbg.open(dbg.verb_modes['verbose'],"debug.log")
    dbg.debug("Debug file\n\n")
    
    sep = '\n----------------------------------------------\n'
    
    print('------ Welcome to the demonstration of -------\n')
    print('               "AFTERSHOQ" \n\n')
    print('       Written by Martin Franckie 2018.')
    print('       Please give credit where credit')
    print('                   is due\n')
    print(sep)
    print('Creating semiconductor materials and alloys:\n')
    print('       ( All CBO relative to GaAs )\n')
    
    # create materials GaAs, AlAs, InAs, alloys:
    
    gaas = GaAs()
    alas = AlAs()
    inas = InAs()
    algaas = AlGaAs(x=0.25)
    ingaasLM = InGaAs()
    inalasLM = AlInAs()
    alingaas = AlInGaAs(x = 0.20)
    
    
    print(str(gaas) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(gaas.params[mp.valdict[val]]))
    
    
    print("\n" + str(alas) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(alas.params[mp.valdict[val]]))
        
    
    print("\n" + str(inas) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(inas.params[mp.valdict[val]]))
    
    print("\n" + str(ingaasLM) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(ingaasLM.params[mp.valdict[val]]))
    
    print("\n" + str(inalasLM) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(inalasLM.params[mp.valdict[val]]))
    
   
    print("\n" + str(algaas) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(algaas.params[mp.valdict[val]]))
        
    print("\n" + str(alingaas) + ":\n")
    for val in mp.valdict:
        print(val + " = " + str(alingaas.params[mp.valdict[val]]))
    
    print(sep)
    print('Creating a structure from generated materials:\n')
    
    # creating a quantum-cascade structure:
    s = Structure()
    
    # set interface roughness parameters for all interfaces:
    s.setIFR( eta = 0.1, lam = 2.0 )
    
    # Add layers:
    s.addLayerWM(8.4,gaas)
    s.addLayerWM(3.1,algaas)
    s.addLayerWM(18.0,gaas)  # <--- this is layer 2!
    s.addLayerWM(1.8,algaas)
    s.addLayerWM(8.4,gaas)
    s.addLayerWM(3.1,algaas)
    s.addLayerWM(18.0,gaas)
    s.addLayerWM(1.8,algaas)
    
    # add a doping layer
    s.addDoping(zi = 2, zf = 2.2, 
                density = 2e17, layerindex = 2)
    
    qcl_2well = EV2416()
    
    print(s)
    print(qcl_2well)
    
    print(sep)
    print('Generating N random structures:\n')
    N = eval(input('N = ?\n'))
    
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
        print(str(st))
        
    print(sep)
    
    proceed = eval(input("Proceed to create directory tree?\nYes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
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
    
    print('Creating directory tree in: ' + path + '.\n')
    print('Please change variable "binpath" to match your binaries!\n')
    
    # make the directory:
    su.mkdir(path)
    
    
    # Define the platform, local or cluster (e. g. Euler cluster at ETH Zuerich)
    pltfm = Local()
    #pltfm = Euler(6,"1:00")
    
    # select model for transport and em raditaion simulations:
    
    # NEGF interface:
    #model = Inegf(binpath,pltfm,numpar,GaAs)
    
    # sewself interface:
    material_list = [gaas,algaas]
    model = Isewself(binpath,pltfm,material_list)
    model.writeStructFile(s, path)
    model.writeParameterFile(path)
    model.writeSewselfPar(path)
    
    # define the merit function as max gain/current density, from a dictionary of merit funcs.:
    model.merit = model.merits.get("(max gain)/(current density)")
    
    print(sep + 'Starting simulations in directory tree?')
    proceed = eval(input("Yes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
        dbg.close()
        exit()
        
    # execute the model for simulating the structures:
    model.runStructures(sg.structures, path)
    
    model.waitforproc(1, str(model) + " running....")
    
    print(sep + 'Gathering results to: ' + path + '/results.log')
    
    # gather the results from the simulation of all structures:
    model.gatherResults(sg.structures, path)
    
    print(sep + 'First results acheived. Proceed with optimization of gain?')
    proceed = eval(input("Yes = 1\nExit = 0"))
    if proceed == 0:
        print(sep + "User exit. Good bye!")
        dbg.close()
        exit()
    
    # collect results from trial points
    x0 = []
    [x0.append( sg.hutil.interp_dist_from_coords( c ) ) for c in coords]
    x0 = numpy.array(x0)
    x0.sort()
    x0 = x0.tolist()
    print(x0)
    y0 = []
    xi = 0
    for i in range(0,len(x0)):
        try:
            y0.append( -float(model.getMerit(sg.structures[i],path)) )
        except( ValueError ):
            del x0[xi]
            xi-=1
        xi +=1
    
    print('Array of trial results:\n' + str(y0))

    # create optimization object
    tol, r, itmax, procmax = 1, 1.3, 10, 2 
    opt = Paraopt(tol,r,itmax,procmax,x0,y0)
    
    conv = opt.minimize(model, sg, path)
    
    
    
    print(sep + 'Program completed! Good bye!')
    dbg.close()
    
    
