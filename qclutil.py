'''
Created on 29 Jan 2018

@author: martin
'''

from classes import Material, MaterialPar as Par

class MaterialUtil():
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
        if name is None:
            name = "GaAs"
        paramsGaAs = []
        Par.initList(paramsGaAs)
        paramsGaAs[0:11] = [0.067,0,1.599,22,0,0.0376,12.9,10.89,-7.2,4730,5317,0.04517]
        paramsGaAs[Par.lattconst] = 5.653
        GaAs = Material(name,paramsGaAs)
        return GaAs
    
    @staticmethod
    def createAlAs(name = None):
        if name is None:
            name = "AlAs"
        paramsAlAs = []
        Par.initList(paramsAlAs)
        paramsAlAs[Par.meff] = 0.15
        paramsAlAs[Par.Ec] = 1.05
        paramsAlAs[Par.Eg] = 3.099
        paramsAlAs[Par.ELO] = 0.050
        paramsAlAs[Par.Ep] = 22
        paramsAlAs[Par.eps0] = 10.06
        paramsAlAs[Par.epsinf] = 8.05 
        AlAs = Material(name,paramsAlAs)
        return AlAs
    
    @staticmethod
    def createInAs(name = None):
        if name is None:
            name = "InAs"
        params = []
        Par.initList(params)
        params[Par.meff] = 0.026
        params[Par.Ec] = -0.892
        params[Par.Eg] = 0.417
        params[Par.ELO] = 0.0298 # Ioffe
        params[Par.Ep] = 21.5
        params[Par.eps0] = 15.1 # Ioffe
        params[Par.epsinf] = 12.3
        params[Par.Vdef] = -5.08
        params[Par.vlong] = 3830
        params[Par.lattconst] = 6.0583
        params[Par.massdens] = 5680
        params[Par.molV] = 0.0556
        InAs = Material(name,params)
        return InAs
    
    @staticmethod
    def createZnO(name = None):
        if name is None:
            name = "ZnO"
        params = []
        Par.initList(params)
        params[Par.meff] = 0.22
        params[Par.Ec] = 0
        params[Par.Eg] = 3.4 # Janotti, van der Walle PRB 2007
        params[Par.ELO] = 0.072
        params[Par.Ep] = 21.5
        params[Par.eps0] = 8.49
        params[Par.epsinf] = 3.72
        params[Par.Vdef] = -2.3 # Janotti, van der Walle PRB 2007
        params[Par.vlong] = 6090 # Bateman JAP 33, 1962
        params[Par.massdens] = 5606
        params[Par.molV] = 0.0145
        params[Par.lattconst] = 5.2
        ZnO = Material(name,params)
        return ZnO
    
    @staticmethod
    def createMgO(name = None):
        if name is None:
            name = "MgO"
        params = []
        Par.initList(params)
        params[Par.meff] = 0.37
        params[Par.Ec] = 1.7 # Janotti, van der Walle PRB 2007
        params[Par.Eg] = 6.3 # Janotti, van der Walle PRB 2007
        params[Par.ELO] = 0.089
        params[Par.Ep] = 21.5
        params[Par.eps0] = 9.6
        params[Par.epsinf] = 2.98
        params[Par.Vdef] = -4.3
        params[Par.lattconst] = 4.21
        MgO = Material(name,params)
        return MgO
    
    @classmethod
    def createAlGaAs(cls,x,name=None):
        if name is None:
            name = "AlGaAs"
        GaAs = cls.createGaAs()
        AlAs = cls.createAlAs()
        A = []
        Par.initList(A)
        A[Par.Ec] = -0.185
        A[Par.Eg] = 0.5
        AlGaAs = Material(name,Material.alloy(AlAs, GaAs, A, x))
        return AlGaAs
    
    @staticmethod
    def createAlInAs(x=None,name=None):
        if name is None:
            name = "AlInAs"
        if x is None:
            x = 0.48
        AlAs = MaterialUtil.createAlAs()
        InAs = MaterialUtil.createInAs()
        A = []
        Par.initList(A)
        A[Par.Eg] = 0.70
        A[Par.meff] = 0.049
        A[Par.Ep] = -4.81
        A[Par.Vdef] = -1.4
        A[Par.Ec] = 0.094 # (= Cvbo + Cgap, Vurgaftman)
        AlInAs = Material(name,Material.alloy(AlAs, InAs, A, x))
        #if x == 0.48:
        #    AlInAs.params[Par.Ec] = 0.52
        return AlInAs
    
    @staticmethod
    def createGaInAs(x=None,name=None):
        if name is None:
            name = "GaInAs"
        if x is None:
            x = 0.47
        GaAs = MaterialUtil.createGaAs()
        InAs = MaterialUtil.createInAs()
        A = []
        Par.initList(A)
        A[Par.Eg] = 0.477
        A[Par.meff] = 0.0091
        A[Par.Ep] = -1.48
        A[Par.Vdef] = 2.61
        A[Par.eps0] = -0.67 # Ioffe
        A[Par.ELO] = 0.002 # fit to Ioffe 34 meV lattice matched
        A[Par.Ec] = 0.060 # (= Cvbo + Cgap, Vurgaftman)
        GaInAs = Material(name,Material.alloy(GaAs, InAs, A, x))
        #GaInAs.params[Par.Ec] = 0 # set Ec relative to this material
        return GaInAs
    
    @staticmethod
    def createZnMgO(x,name=None):
        if name is None:
            name = "Zn_"+str(x)+"Mg_"+str(1-x)+"O"
        ZnO = MaterialUtil.createZnO()
        MgO = MaterialUtil.createMgO()
        # Bowing paremeters:Unknown
        A = []
        Par.initList(A)
        ZnMgO = Material(name,Material.alloy(ZnO,MgO,A,x))
        return ZnMgO