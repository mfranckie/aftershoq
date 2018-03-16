'''
Created on 26 Jan 2018

@author: martin
'''

from classes import MaterialPar as Par
from systemutil import SystemUtil as su
import time
from runplatf import Local

class Interface(object):
    
    merits = {"max gain":0, 
              "(max gain)/(current density)":1, 
              "inversion":2, 
              "output power":3, 
              "wall plug efficiency":4,
              "threshold current density":5,
              "estimated gain":6,
              "custom figure of merit":7}
    
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