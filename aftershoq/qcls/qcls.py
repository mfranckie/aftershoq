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
        
        self.addLayerWM(3.1, barrier)
        self.addLayerWM(8.5, well)
        self.addLayerWM(1.8, barrier)
        self.addLayerWM(8.7, well)
        self.addLayerWM(3.0, well) # <----- # 4 doped to 1.5*10^17 cm^-3 (4.5*10^10 cm^-2)
        self.addLayerWM(6.0, well)
        
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
        
        self.addLayerWM(5.5, barrier)
        self.addLayerWM(11.0, well)
        self.addLayerWM(1.8, barrier)
        self.addLayerWM(11.5, well)
        self.addLayerWM(3.8, barrier) 
        self.addLayerWM(9.4, well)
        self.addLayerWM(4.2, barrier) 
        self.addLayerWM(18.4, well) # <----- # 7 doped to 2*10^16 cm^-3
        
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
        
        self.addLayerWM(4.3, barrier)
        self.addLayerWM(8.9, well)
        self.addLayerWM(2.46, barrier)
        self.addLayerWM(8.15, well)
        self.addLayerWM(4.1, barrier) 
        self.addLayerWM(5.5, well) 
        self.addLayerWM(5.0, well) # <----- # 6 doped to 6*10^16 cm^-3 (4.5*10^10 cm^-2)
        self.addLayerWM(5.5, well)
        
        idop= 6
        vdop = 6e16 # cm^-3
        
        self.addDoping(0, self.layers[idop].width, vdop, idop)
    

class N471(Structure):
    '''THz QCL emitting at 3.74 THz "single quantum well active region".
    Published in: Scalari et al., Appl. Phys. Lett. 91, 032103 (2007)
    '''
    def __init__(self):
        Structure.__init__(self)
        
        self.setIFR(0.1, 10)
        well = GaAs()
        barrier = AlGaAs(x = 0.15)
        
        self.addLayerWM(4.7, barrier)
        self.addLayerWM(28, well)
        self.addLayerWM(2.3, barrier)
        self.addLayerWM(18, well)
        self.addLayerWM(2.3, barrier) 
        self.addLayerWM(16.5, well)
        self.addLayerWM(2.3, barrier) 
        self.addLayerWM(16.0, well) # <----- # 7 doped to 2.4*10^16 cm^-3 (3.84*10^10 cm^-2)
        self.addLayerWM(2.3, barrier) 
        self.addLayerWM(15.5, well)
        
        idop= 7
        vdop = 2.4e16 # cm^-3
        
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
        
        
        self.addLayerWM(3.1, alinas_s)
        self.addLayerWM(2.5, gainas_s)
        self.addLayerWM(0.6, alinas_s)
        self.addLayerWM(5.7, gainas_s)
        self.addLayerWM(0.7, alinas_s)
        self.addLayerWM(5.5, gainas_s)
        self.addLayerWM(1.2, alinas_s)
        self.addLayerWM(4.6, gainas_s)
        self.addLayerWM(1.1, alinas_s)
        self.addLayerWM(4.5, gainas_s) # <----- #9  Doped to 0.120276*10^18 cm^-3
        self.addLayerWM(1.4, alinas_s) # <----- #10 Doped to 0.120276*10^18 cm^-3 
        self.addLayerWM(4.0, gainas_s) # <----- #11 Doped to 0.120276*10^18 cm^-3
        self.addLayerWM(1.5, alinas_s)
        self.addLayerWM(3.4, gainas_s)
        self.addLayerWM(1.7, alinas_s)
        self.addLayerWM(3.5, gainas_s)
        
        dop = 0.10101e18
        
        self.addDoping(0, 4.5, dop, 9)
        self.addDoping(0, 1.4, dop, 10)
        self.addDoping(0, 4.0, dop, 11)
        
        
class EV2016(Structure):
    '''
    Strained 8.5 micron design, seed for the genetically optimized design EV2017.
    Published in thesis of J. Wolf (ETH Zuerich, 2017)
    and
    Optics Express vol. 20(22) p.24272 (2012)
    '''
        
    def __init__(self):
        Structure.__init__(self)
        
        self.setIFR(0.1, 10)
        
        alinas_s = AlInAs(x = 0.64)
        # to account for strain (see thesis Tobias Gresch ETH Zürich 2009)
        alinas_s.params[mp.Ec] = 0.72
        gainas_s = InGaAs(x = 0.58)
        gainas_s.params[mp.Ec] = 0.0
        
        
        self.addLayerWM(3.1, alinas_s)
        self.addLayerWM(2.5, gainas_s)
        self.addLayerWM(0.6, alinas_s)
        self.addLayerWM(5.7, gainas_s)
        self.addLayerWM(0.7, alinas_s)
        self.addLayerWM(5.5, gainas_s)
        self.addLayerWM(1.2, alinas_s)
        self.addLayerWM(4.6, gainas_s)
        self.addLayerWM(1.1, alinas_s)
        self.addLayerWM(4.5, gainas_s) # <----- #9  Doped to 0.101*10^18 cm^-3
        self.addLayerWM(1.4, alinas_s) # <----- #10 Doped to 0.101*10^18 cm^-3 
        self.addLayerWM(4.0, gainas_s) # <----- #11 Doped to 0.101*10^18 cm^-3
        self.addLayerWM(1.5, alinas_s)
        self.addLayerWM(3.4, gainas_s)
        self.addLayerWM(1.7, alinas_s)
        self.addLayerWM(3.5, gainas_s)
        
        dop = 0.120276e18
        
        self.addDoping(0, 4.5, dop, 9)
        self.addDoping(0, 1.4, dop, 10)
        self.addDoping(0, 4.0, dop, 11)
        
class EV2017(Structure):
    '''
    Genetically optimized 8.5 micron design. 
    Published in thesis of J. Wolf (ETH Zuerich, 2017)
    '''
    
    def __init__(self):
        Structure.__init__(self)
        
        self.setIFR(0.1, 10)
        
        alinas_s = AlInAs(x = 0.64)
        # to account for strain (see thesis Tobias Gresch ETH Zürich 2009)
        alinas_s.params[mp.Ec] = 0.72
        gainas_s = InGaAs(x = 0.58)
        gainas_s.params[mp.Ec] = 0.0
        
        
        self.addLayerWM(3.1, alinas_s)
        self.addLayerWM(2.52, gainas_s)
        self.addLayerWM(1.23, alinas_s)
        self.addLayerWM(5.77, gainas_s)
        self.addLayerWM(0.74, alinas_s)
        self.addLayerWM(5.0, gainas_s)
        self.addLayerWM(1.01, alinas_s)
        self.addLayerWM(4.49, gainas_s)
        self.addLayerWM(1.27, alinas_s)
        self.addLayerWM(3.79, gainas_s) # <----- #9  Doped to 0.120276*10^18 cm^-3
        self.addLayerWM(1.29, alinas_s) # <----- #10 Doped to 0.120276*10^18 cm^-3 
        self.addLayerWM(3.23, gainas_s) # <----- #11 Doped to 0.120276*10^18 cm^-3
        self.addLayerWM(1.60, alinas_s)
        self.addLayerWM(2.89, gainas_s)
        self.addLayerWM(1.89, alinas_s)
        self.addLayerWM(3.01, gainas_s)
        
        dop = 0.120276e18
        
        self.addDoping(0, 3.79, dop, 9)
        self.addDoping(0, 1.29, dop, 10)
        self.addDoping(0, 3.23, dop, 11)
        
        