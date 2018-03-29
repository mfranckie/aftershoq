'''
Created on 29 Jan 2018

@author: martin
'''

from structure.materials import AlGaAs,GaAs,AlAs,ZnO,MgO,InAs,InGaAs,InAlAs,ZnMgO

class MaterialUtil:
    '''
    Help function for creating predefined materials
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.GaAs = self.createGaAs()
        self.AlAs = self.createAlAs()
        self.InAs = self.createInAs()
        self.GaInAs = self.createGaInAs()
        self.AlInAs = self.createAlInAs()
    
    @staticmethod
    def createGaAs(name=None):
        return GaAs(name)
    
    @staticmethod
    def createAlAs(name = None):
        return AlAs(name)
    
    @staticmethod
    def createInAs(name = None):
        return InAs(name)
    
    @staticmethod
    def createZnO(name = None):
        return ZnO(name)
    
    @staticmethod
    def createMgO(name = None):
        return MgO(name)
    
    @classmethod
    def createAlGaAs(cls,x,name=None):
        return AlGaAs(name,x)
    
    @staticmethod
    def createAlInAs(x=None,name=None):
        return InAlAs(name,x)
    
    @staticmethod
    def createGaInAs(x=None,name=None):
        return InGaAs(name,x)
    
    @staticmethod
    def createZnMgO(x,name=None):
        return ZnMgO(name,x)