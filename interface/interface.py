'''
Created on 26 Jan 2018

@author: martin
'''

class Interface(object):
    
    merits = {"max gain":0, 
              "(max gain)/(current density)":1, 
              "inversion":2, 
              "output power":3, 
              "wall plug efficiency":4,
              "threshold current density":5,
              "estimated gain":6,
              "custom figure of merit":7
              }
    
    numpar = {"efield0":0,
              "defield":0.001,
              "Nefield":1,
              "omega0":0,
              "domega":0.001,
              "Nomega":21,
              "efac0": 0,
              "defac": 0.010,
              "Nefac": 1,
              "Tlattice": 77,
              "Te": 125,
              "Nstates":5,
              "Nper":1,
              "maxits":50,
              "NE":1000,
              "Nk":800,
              "Nz":400,
              "Nq":400,
              "use-ifr":True,
              "use-alloy":True,
              "use-acc": True,
              "use-LO": True,
              "use-TO": True,
              "use-imp": True,
              "use-e-e": False,
              "use-poisson": True
              }
    
    def __init__(self,binpath,pltfm):
        self.binpath = binpath
        self.pltfm = pltfm
        self.merit = Interface.merits.get("max gain")

    def runStructures(self,structures,numpar,path):
        pass
    
    def gatherResults(self,structures,path):
        pass
    
    def waitforproc(self, delay, message = None):
        pass