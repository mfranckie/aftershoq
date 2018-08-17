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
    Amanti et al., New J. Phys. 11, 125022 (2009)
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
        self.addLayerMW(18.4, well) # <----- # 7 doped to 2*10^16 cm^-3
        
        idop= 7
        vdop = 2e16 # cm^-3
        
        self.addDoping(0, self.layers[idop].width, vdop, idop)
        
class Fathololoumi2012(Structure):
    '''Record THz QCL from Fathololoumi et al., Optics Express 20, 3866 (2012)
    '''
    def __init__(self):
        Structure.__init__(self)
        
        self.setIFR(0.1, 10)
        well = GaAs()
        barrier = AlGaAs(x = 0.15)
        
        self.addLayerMW(4.3, barrier)
        self.addLayerMW(8.9, well)
        self.addLayerMW(2.46, barrier)
        self.addLayerMW(8.15, well)
        self.addLayerMW(4.1, barrier) 
        self.addLayerMW(5.5, well) 
        self.addLayerMW(5.0, well) # <----- # 6 doped to 6*10^16 cm^-3 (4.5*10^10 cm^-2)
        self.addLayerMW(5.5, well)
        
        idop= 6
        vdop = 6e16 # cm^-3
        
        self.addDoping(0, self.layers[idop].width, vdop, idop)
    
        
class EV1907(Structure):
    '''
    Strained 8.5 micron design. Published in thesis of J. Wolf (ETH Zuerich, 2017)
    '''
    
    def __init__(self):
        Structure.__init__(self)
        
        self.setIFR(0.1, 10)
        
        alinas_s = AlInAs(x = 0.64)
        gainas_s = InGaAs(x = 0.58)
        
        
        self.addLayerMW(3.1, alinas_s)
        self.addLayerMW(2.5, gainas_s)
        self.addLayerMW(0.6, alinas_s)
        self.addLayerMW(5.7, gainas_s)
        self.addLayerMW(0.7, alinas_s)
        self.addLayerMW(5.5, gainas_s)
        self.addLayerMW(1.2, alinas_s)
        self.addLayerMW(4.6, gainas_s)
        self.addLayerMW(1.1, alinas_s)
        self.addLayerMW(4.5, gainas_s) # <----- #9  Doped to 0.120276*10^18 cm^-3
        self.addLayerMW(1.4, alinas_s) # <----- #10 Doped to 0.120276*10^18 cm^-3 
        self.addLayerMW(4.0, gainas_s) # <----- #11 Doped to 0.120276*10^18 cm^-3
        self.addLayerMW(1.5, alinas_s)
        self.addLayerMW(3.4, gainas_s)
        self.addLayerMW(1.7, alinas_s)
        self.addLayerMW(3.5, gainas_s)
        
        dop = 0.10101e18
        
        self.addDoping(0, 4.5, dop, 9)
        self.addDoping(0, 1.4, dop, 10)
        self.addDoping(0, 4.0, dop, 11)
        
class EV2017(Structure):
    '''
    Genetically optimized 8.5 micron design. Published in thesis of J. Wolf (ETH Zuerich, 2017)
    '''
    
    def __init__(self):
        Structure.__init__(self)
        
        self.setIFR(0.1, 10)
        
        alinas_s = AlInAs(x = 0.64)
        gainas_s = InGaAs(x = 0.58)
        
        
        self.addLayerMW(3.1, alinas_s)
        self.addLayerMW(2.52, gainas_s)
        self.addLayerMW(1.23, alinas_s)
        self.addLayerMW(5.77, gainas_s)
        self.addLayerMW(0.74, alinas_s)
        self.addLayerMW(5.0, gainas_s)
        self.addLayerMW(1.01, alinas_s)
        self.addLayerMW(4.49, gainas_s)
        self.addLayerMW(1.27, alinas_s)
        self.addLayerMW(3.79, gainas_s) # <----- #9  Doped to 0.120276*10^18 cm^-3
        self.addLayerMW(1.29, alinas_s) # <----- #10 Doped to 0.120276*10^18 cm^-3 
        self.addLayerMW(3.23, gainas_s) # <----- #11 Doped to 0.120276*10^18 cm^-3
        self.addLayerMW(1.60, alinas_s)
        self.addLayerMW(2.89, gainas_s)
        self.addLayerMW(1.89, alinas_s)
        self.addLayerMW(3.01, gainas_s)
        
        dop = 0.120276e18
        
        self.addDoping(0, 3.79, dop, 9)
        self.addDoping(0, 1.29, dop, 10)
        self.addDoping(0, 3.23, dop, 11)
        
        