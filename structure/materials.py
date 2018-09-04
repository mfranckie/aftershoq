'''
Created on 29 Mar 2018

@author: Martin Franckie

This module contains pre-defined semiconductor materials and
binary, ternary, as well as quaternary alloys.

For materials where the alloy composition does not imply a
simple bowing relationship for some parameters, the parent
updateAlloy() and copy() functions are overridden.

Note: parameters which are not explicitly altered are set to 0 by default!

Naming convention: For group XX-YY materials, the XX material comes before 
the YY material. For alloys XX_xYY_(1-x), the XX material comes before the
YY material.

References:
[Vurgaftman2001] Vurgaftman et al., Appl. Phys. Rev. 89, 5815-5875 (2001)
[Ioffe] http://matprop.ru and sources therein.
[Janotti] Janotti and van der Walle, Phys. Rev. B 75, 121201 (2007)
[Bateman] Bateman et al., J. Appl. Phys. 33, 3309 (1962)
'''

from structure.classes import Material
import structure.matpar as mp

# Binaries:

class GaAs(Material):
    '''GaAs, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe].
    '''
    
    def __init__(self,name = None):
        if name is None:
            name = "GaAs"
        paramsGaAs = []
        mp.initList(paramsGaAs)
        paramsGaAs[mp.meff] = 0.067
        paramsGaAs[mp.Ec] = 0.0
        paramsGaAs[mp.Eg] = 1.519
        paramsGaAs[mp.ELO] = 0.0376
        paramsGaAs[mp.Ep] = 22
        paramsGaAs[mp.eps0] = 12.9
        paramsGaAs[mp.epsinf] = 10.89
        paramsGaAs[mp.Vdef] = -7.17
        paramsGaAs[mp.vlong] = 4730
        paramsGaAs[mp.massdens] = 5317
        paramsGaAs[mp.molV] = 0.04517
        paramsGaAs[mp.lattconst] = 5.653
        super(GaAs,self).__init__(name,paramsGaAs)
        
    def copy(self):
        return GaAs(self.name)

class AlAs(Material):
    '''AlAs, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe].
    '''
    def __init__(self,name = None):
        if name is None:
            name = "AlAs"
        paramsAlAs = []
        mp.initList(paramsAlAs)
        paramsAlAs[mp.meff] = 0.15
        paramsAlAs[mp.Ec] = 0.99
        paramsAlAs[mp.Eg] = 3.099
        paramsAlAs[mp.ELO] = 0.050
        paramsAlAs[mp.Ep] = 22
        paramsAlAs[mp.eps0] = 10.06
        paramsAlAs[mp.epsinf] = 8.05 
        super(AlAs,self).__init__(name,paramsAlAs)
        
    def copy(self):
        return AlAs(self.name)

class InAs(Material):
    '''InAs, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe].
    '''    
    def __init__(self,name = None):
        
        if name is None:
            name = "InAs"
            params = []
            mp.initList(params)
            params[mp.meff] = 0.026
            params[mp.Ec] = -0.892
            params[mp.Eg] = 0.417
            params[mp.ELO] = 0.0298 # Ioffe
            params[mp.Ep] = 21.5
            params[mp.eps0] = 15.1 # Ioffe
            params[mp.epsinf] = 12.3
            params[mp.Vdef] = -5.08
            params[mp.vlong] = 3830
            params[mp.lattconst] = 6.0583
            params[mp.massdens] = 5680
            params[mp.molV] = 0.0556
            super(InAs,self).__init__(name,params)
            
    def copy(self):
        return InAs(self.name)

class ZnO(Material):
    '''ZnO, parameters from [Janotti2007] and [Bateman1962].
    Note that at this time, the state of knowledge about this
    material is very limited.'''
    def __init__(self,name = None):
        if name is None:
            name = "ZnO"
        params = []
        mp.initList(params)
        params[mp.meff] = 0.22
        params[mp.Ec] = 0
        params[mp.Eg] = 3.4 # [Janotti2007]
        params[mp.ELO] = 0.072
        params[mp.Ep] = 21.5
        params[mp.eps0] = 8.49
        params[mp.epsinf] = 3.72
        params[mp.Vdef] = -2.3 # [Janotti2007]
        params[mp.vlong] = 6090 # Bateman JAP 33, 1962
        params[mp.massdens] = 5606
        params[mp.molV] = 0.0145
        params[mp.lattconst] = 5.2
        super(ZnO,self).__init__(name,params)
        
    def copy(self):
        return ZnO(self.name)
        
class MgO(Material):
    '''MgO, parameters from [Janotti2007].
    Note that at this time, the state of knowledge about this
    material is very limited.'''
    def __init__(self,name = None):
        
        if name is None:
            name = "MgO"
        params = []
        mp.initList(params)
        params[mp.meff] = 0.37
        params[mp.Ec] = 1.7 # Janotti, van der Walle PRB 2007
        params[mp.Eg] = 6.3 # Janotti, van der Walle PRB 2007
        params[mp.ELO] = 0.089
        params[mp.Ep] = 21.5
        params[mp.eps0] = 9.6
        params[mp.epsinf] = 2.98
        params[mp.Vdef] = -4.3
        params[mp.lattconst] = 4.21
        super(MgO,self).__init__(name,params)
        
    def copy(self):
        return MgO(self.name)
    
class MgOzoterac(Material):
    '''MgO, parameters used in the Zoterac (ERC) project.'''
    
    def __init__(self,name = 'MgOzoterac'):
        params = []
        mp.initList(params)
        params[mp.meff] = 0.22
        params[mp.Ec] = 1.3 # Deliverable D1.3 1/9/2017
        params[mp.Eg] = 5.3 # Deliverable D1.3 1/9/2017
        params[mp.ELO] = 0.089
        params[mp.Ep] = 21.5
        params[mp.eps0] = 9.6
        params[mp.epsinf] = 2.98
        params[mp.Vdef] = -4.3
        params[mp.lattconst] = 4.21
        super(MgOzoterac,self).__init__(name,params)

    def copy(self):
        return MgOzoterac(self.name)
        
# Ternaries:

class AlGaAs(Material):
    '''Al_xGa_1-xAs. Bowing parameters from [Vurgaftman2001].'''

    def __init__(self, name = None, x=None):
        if name is None:
            name = "AlGaAs"
        mat1 = AlAs()
        mat2 = GaAs()
        C = []
        mp.initList(C)
        C[mp.Eg] = -0.127 + 1.310*x
        super(AlGaAs,self).__init__(name, [], mat1, mat2, C, x)
            
    def updateAlloy(self,x):
        super(AlGaAs,self).updateAlloy(x)
        if x < 0.42:
            self.params[mp.Ec] = 0.831*x
        else:
            self.params[mp.Ec] = 0.332 + 0.054*x
            
    def copy(self):
        return AlGaAs(self.name,self.x)
    
    
    
class InGaAs(Material):
    '''In_xGa_1-xAs. Bowing parameters from [Vurgaftman2001] and [Ioffe].'''
    def __init__(self,name = None, x = None):
        if name is None:
            name = "GaInAs"
        if x is None:
            x = 0.47
        mat1 = GaAs()
        mat2 = InAs()
        A = []
        mp.initList(A)
        A[mp.Eg] = 0.477
        A[mp.meff] = 0.0091
        A[mp.Ep] = -1.48
        A[mp.Vdef] = 2.61
        A[mp.eps0] = -0.67 # Ioffe
        A[mp.ELO] = 0.002 # fit to Ioffe 34 meV lattice matched
        A[mp.Ec] = 0.060 # (= Cvbo + Cgap, Vurgaftman)
        super(InGaAs,self).__init__(name,[],mat1, mat2, A, x)
        
        if x == 0.47:
            self.params[mp.Eg] = 0.789
        
    def copy(self):
        return InGaAs(self.name,self.x)
        
class AlInAs(Material):
    '''Al_xIn_1-xAs. Bowing parameters from [Vurgaftman2001].'''
    def __init__(self,name = None, x = None):
        
        if name is None:
            name = "AlInAs"
        
        mat1 = AlAs()
        mat2 = InAs()
        A = []
        mp.initList(A)
        A[mp.Eg] = 0.70
        A[mp.meff] = 0.049
        A[mp.Ep] = -4.81
        A[mp.Vdef] = -1.4
        A[mp.Ec] = 0.094 # (= Cvbo + Cgap, Vurgaftman)
        if x is None:
            x = 0.48
            
        super(AlInAs,self).__init__(name,[],mat1, mat2, A, x)
        
        if x == 0.48:
            self.params[mp.Ec] = 0.032294
            self.params[mp.meff] = 0.080
            self.params[mp.Eg] = 1.404
        
    def copy(self):
        return AlInAs(self.name,self.x)
        


class ZnMgO(Material):
    '''Zn_xMg_1-xO. Only linear mixing.'''
    def __init__(self,name = None,x = None):
        if name is None:
            name = "Zn_"+str(x)+"Mg_"+str(1-x)+"O"
        mat1 = ZnO()
        mat2 = MgO()
        # Bowing paremeters:Unknown
        A = []
        mp.initList(A)
        super(ZnMgO,self).__init__(name,[],mat1,mat2,A,x)
        
    def copy(self):
        return ZnMgO(self.name,self.x)
    
    
# Quqternaries;
    
class AlInGaAs(Material):
        '''Quaternary alloy composed of AlInAs and GaInAs
        Default: Lattice matched to InP (Al_0.48InAs)_x (In_0.53GaAs)_(1-x)
        Material parameters from Ohtani APL (2013)
        Note: General case is not implemented yet!
        
        Using a linear interpolation of the alloy scattering potential.
        ''' 
        
        def __init__(self, x, name = None, y = None, z = None):
            
            if name is None:
                name = "Al_"+str(x) + "Ga_"+ str(1-x) + "InAs"
            
            if y is None:
                y = 0.48
            elif y != 0.48:
                print "General case is not implemented, use y = 0.48 or None!"
                
            if z is None:
                z = 0.53
            elif z != 0.53:
                print "General case is not implemented, use z = 0.53 or None!"
            
            mat1 = AlInAs()
            mat2 = InGaAs()
            
            A = []
            mp.initList(A)
            A[mp.Eg] = 0.22
            A[mp.meff] = -0.016
            super(AlInGaAs,self).__init__(name,[],mat1,mat2,A,x)
            
        def copy(self):
            return AlInGaAs(self.x, self.name)
        
        def updateAlloy(self, x):
            Material.updateAlloy(self, x)
            self.params[mp.Valloy] = self.mat1.params[mp.Valloy]*x + \
                self.mat2.params[mp.Valloy]*(1-x)
            self.params[mp.Ec] = 0.73*(0.712*x - 0.22*x*(1-x)) + self.mat2.params[mp.Ec]
        