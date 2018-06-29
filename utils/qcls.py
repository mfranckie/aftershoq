'''
Created on 19 Jun 2018

@author: martin

This file contains public structures for reference. Please do not add confidential
layer sequences to this list! Please give a cited reference to each structure, 
where possible.
'''

from structure.classes import Structure
from structure.materials import *

class EV2416(Structure):
    '''
    Record (192 K) 2-well THz QCL
    Franckie et al. Appl. Phys. Lett. 112, 021104 (2018)
    '''
    
    def __init__(self):
        Structure.__init__(self)
        
        self.setIFR(0.1, 10)
        well = GaAs()
        barrier = AlGaAs(x = 0.25)
        
        self.addLayerMW(3.1, barrier)
        self.addLayerMW(8.5, well)
        self.addLayerMW(1.8, barrier)
        self.addLayerMW(8.7, well)
        self.addLayerMW(3.0, well) # <----- # 4 doped to 1.5*10^17 cm^-3 (4.5*10^10 cm^-2)
        self.addLayerMW(6.0, well)
        
        idop= 4
        vdop = 1.5e17 # cm^-3
        
        self.addDoping(0, self.layers[idop].width, vdop, idop)
        
class EV1157(Structure):
    '''
    Amanit et al. New J. Phys. 11, 125022 (2009)
    '''
    
    def __init__(self):
        Structure.__init__(self)
        
        self.setIFR(0.1, 10)
        well = GaAs()
        barrier = AlGaAs(x = 0.15)
        
        self.addLayerMW(5.5, barrier)
        self.addLayerMW(11.0, well)
        self.addLayerMW(1.8, barrier)
        self.addLayerMW(11.5, well)
        self.addLayerMW(3.8, barrier) 
        self.addLayerMW(9.4, well)
        self.addLayerMW(4.2, barrier) 
        self.addLayerMW(18.4, well) # <----- # 7 doped to 1.5*10^17 cm^-3 (4.5*10^10 cm^-2)
        
        idop= 7
        vdop = 2e16 # cm^-3
        
        self.addDoping(0, self.layers[idop].width, vdop, idop)