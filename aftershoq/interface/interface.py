'''
Created on 26 Jan 2018

@author: Martin Franckie

Module containing Interface class. This class is intended as parent class to
model interfaces, which should implement the following methods:

runStructures( structures, path )   :  Run structures in list "structures".
                                       "Path" is the working directory which
                                       contains the tree of structures.
gatherResults( structures, path )   :  Prepare and gather model results, 
                                       put the results in the respective 
                                       structure as new attributes.
waitforproc()                       :  Wait for all processes to finish.
__init__()                          :  Optionally add model-specific arguments

If the model uses additional numerical parameters, define one or more new 
dictionaries and in __init__() call self.numpar.update( _new_dict_ )

'''

import time
from aftershoq.structure import Structure

class Interface(object):
    
    merits = {"max gain":0, 
              "(max gain)/(current density)":1, 
              "inversion":2, 
              "output power":3, 
              "wall plug efficiency":4,
              "threshold current density":5,
              "estimated gain":6,
              "Chi2":7,
              "custom figure of merit":8
              }
    
    # Energies in meV, temperatures in K
    numpar = {"efield0":0,
              "defield":0.001,
              "Nefield":1,
              "omega0":0.001,
              "domega":0.001,
              "Nomega":1,
              "efac0": 0,
              "defac": 0.010,
              "Nefac": 1,
              "Tlattice": 77,
              "Te": 125,
              "Nstates":5,
              "Nper":1,
              "maxits":50,
              "NE":1000,
              "Emax":1.0,
              "Nk":800,
              "Ekmax":1.0,
              "Nz":400,
              "Nq":400,
              "use-ifr":True,
              "use-alloy":True,
              "use-acc": True,
              "use-LO": True,
              "use-TO": True,
              "use-imp": True,
              "use-e-e": False,
              "use-poisson": True,
              "use-poisson": True
              }
    
    def __init__(self,binpath,pltfm):
        '''Constructor.
        binpath : path to model binary files
        pltfm : Platf object
        
        '''
        self.binpath = binpath
        self.pltfm = pltfm
        self.merit = Interface.merits.get("max gain")

    def runStructures(self,structures,path,runprog=True):
        '''Run simulations for all structures in the given structure list with
        the base path "path". This method dispatches all processes and returns
        the user has to wait for processes to finish before accessing results.
        
        Stores started processes in self.processes
        '''
        pass
    
    def gatherResults(self,structures,path,pathresults=None, runprog=None):
        '''Write results to pathresults/results.log and run hdiag and bandplot
        in pathwd/s.dirname/self.datpath/eins/x/ for each i and x. Stores WS 
        resutls as a new attribute levels[directory][WS level][data field] in 
        each Structure object in the list structures.
        '''
        pass
    
    def waitforproc(self, delay, message = None):
        '''Blocks execution until all processes in self.processes are 
        finished.
        '''
        
        pactive = True
        while pactive:
            if message is not None:
                print(message)
            pactive = False
            for p in self.processes:
                if self.pltfm.jobstatus(p):
                    pactive=True
                    #break
            time.sleep(delay)
        # close all finished processes:
        for p in self.processes:
            p.wait()
            del p
    
    def getMerit(self,structure,path):
        '''Returns the merit function evaluated for the Structure structure,
        with base path "path". 
        '''
        pass
    
    @classmethod
    def loadStructure(cls, resultpath, origs, sid):
        """
        Loads and returns the structure with id "sid" based on the "origs" Structure,
        from the results file at "resultpath/results.log".
        Returns None if no structure with the given id was found.
        """
        
        s = Structure(origs)
                      
        Nl = len(origs.layers)
        Ndop = len(origs.dopings)
        
        found = False
        
        with open(resultpath + "/results.log") as resfile:
            try:
                for line in resfile:
                    line = line.split()
                    if line[0] == str(sid):
                        found=True
                        break
            except(IndexError):
                pass
            
            if not found:
                return None
            
            for il in range(Nl):
                w = float( line[il+2] )
                x = line[il+2+Ndop*3+Nl]
                try:
                      x = float(x)
                      s.layers[il].material.updateAlloy(x)
                except(ValueError):
                    pass # x = 'None'
                s.layers[il].width = w
                
            for idop in range(Ndop):
                for i in range(3):
                    s.dopings[idop][i] = float(line[Nl + 2 + idop*3 + i])
                
        return s
    
    @classmethod
    def loadAllStructures(cls, resultpath, origs):



        Nl = len(origs.layers)
        Ndop = len(origs.dopings)
        Structure.sid = 0
        structures = []

        with open(resultpath + "/results.log") as resfile:
            for line in resfile:
                
                line = line.split()

                try:
                    sid = int(line[0])
                except(ValueError, IndexError):
                    continue # comment or empty line, skip
                s = Structure(origs)
                for il in range(Nl):
                    w = float( line[il+2] )
                    x = line[il+2+Ndop*3+Nl]
                    try:
                          x = float(x)
                          s.layers[il].material.updateAlloy(x)
                    except(ValueError):
                        pass # x = 'None'
                    s.layers[il].width = w

                for idop in range(Ndop):
                    for i in range(3):
                        s.dopings[idop][i] = float(line[Nl + 2 + idop*3 + i])

                structures.append(s)

        return structures
                
        