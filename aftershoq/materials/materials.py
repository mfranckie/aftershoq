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

from aftershoq.structure import Material

class Si(Material):
    '''
    Ref: K. Driscoll and R. Paiella, J. Appl. Phys. 102, 093103 (2007)
    '''
    
    def __init__(self,name = None):
        if name is None:
            name = "Si"
        paramsSi = Material.params_dict.copy()
        paramsSi["meff"] = 0 # check
        paramsSi["Ec"] = 0.0 # check
        paramsSi["Eg"] = 3.37 # <--- Probably assumes high Ge content
        paramsSi["EDel"] = 1.155 # <--- Probably assumes high Ge content
        paramsSi["EL"] = 2.01 # <--- Probably assumes high Ge content
        paramsSi["ELO"] = 0 # check
        paramsSi["Ep"] =0 # check
        paramsSi["eps0"] = 0 # check
        paramsSi["epsinf"] = 0 # check
        paramsSi["ac"] = 0 # Only know difference
        paramsSi["av"] = 0 #  Only know difference
        paramsSi["c11"] = 16.75 # (10^6 N cm^{-2})
        paramsSi["c12"] = 6.5
        paramsSi["c44"] = 0 # check
        paramsSi["vlong"] = 0 # check
        paramsSi["massdens"] = 0 # check
        paramsSi["molV"] = 0 # check
        paramsSi["lattconst"] = 5.431
        super(Si,self).__init__(name,paramsSi)
        
    def copy(self):
        return Si(self.name)
    
class Ge(Material):
    '''
    Ref: K. Driscoll and R. Paiella, J. Appl. Phys. 102, 093103 (2007)
    '''
    
    def __init__(self,name = None):
        if name is None:
            name = "Ge"
        paramsGe = Material.params_dict.copy()
        paramsGe["meff"] = 0.12
        paramsGe["Ec"] = 0.0 # check
        paramsGe["Eg"] = 0.89 # <--- Probably assumes high Ge content
        paramsGe["EDel"] = 0.931 # <--- Probably assumes high Ge content
        paramsGe["EL"] = 0.74 # <--- Probably assumes high Ge content
        paramsGe["ELO"] = 0 # Check
        paramsGe["Ep"] = 0 # check
        paramsGe["eps0"] = 0 # check
        paramsGe["epsinf"] = 0 # check
        paramsGe["ac"] = -7.17 # only know difference
        paramsGe["av"] = 1.16 # check
        paramsGe["c11"] = 13.15
        paramsGe["c12"] = 4.94
        paramsGe["c44"] = 0 # check
        paramsGe["vlong"] = 0 # check
        paramsGe["massdens"] = 0 # check
        paramsGe["molV"] = 0 # check
        paramsGe["lattconst"] = 5.657
        super(Ge,self).__init__(name,paramsGe)
        
    def copy(self):
        return Ge(self.name)

# Binaries:
        
class SiGe(Material):
    '''
    Ref: K. Driscoll and R. Paiella, J. Appl. Phys. 102, 093103 (2007)
    Ref for Gamma: Phys. Rev. B 82, 205317 (2010)
    
    Exception with definition!!
    Si_(1-x)Ge_x
    '''

    def __init__(self, name = None, x=0.):
        if name is None:
            name = "Si_" + str(1-x) + "Ge_" + str(x)
        mat1 = Si()
        mat2 = Ge()
        C = Material.params_dict.copy()
        # Lattice constant from DJ Paul
        C["lattconst"] = 5.431 + 0.1992*x + 0.02733*x**2
        # Unstrained bandgaps
        C["Eg"] = 3.37 - 2.48*x
        C["EL"] = 2.01 - 1.27*x
        C["EDel"] = 1.155-0.43*x+0.206*x**2
        super(SiGe,self).__init__(name, Material.params_dict.copy(), 
                                    mat1, mat2, x)
            
#    def updateAlloy(self,x,reset_strain = False):
#        self.C["Eg"] = 3.37 - 2.48*x
#        self.C["EL"] = 2.01 - 1.27*x
#        self.C["EDel"] = 1.155-0.43*x+0.206*x**2
#        super(SiGe,self).updateAlloy(x,reset_strain = reset_strain)
#        if self.CBO_Yi:
#            if x < 0.42:
#                self.params["Ec"] = 0.831*x
#            else:
#                #self.params["Ec"] = 0.332 + 0.054*x #indirect gap
#                self.params["Ec"] = -0.115 + 1.105*x  #direct gap
            
    def copy(self):
        return SiGe(self.name,self.x)

class GaAs(Material):
    '''GaAs, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe]. Elastic constants and hydrostatic potentials from 
    [vandeWalle1989].
    '''
    
    def __init__(self,name = None):
        if name is None:
            name = "GaAs"
        paramsGaAs = Material.params_dict.copy()
        paramsGaAs["meff"] = 0.067
        paramsGaAs["Ec"] = 0.0
        paramsGaAs["Eg"] = 1.519
        paramsGaAs["EX"] = 1.981
        paramsGaAs["ELO"] = 0.0376
        paramsGaAs["Ep"] = 22
        paramsGaAs["eps0"] = 12.9
        paramsGaAs["epsinf"] = 10.89
        paramsGaAs["ac"] = -7.17
        paramsGaAs["av"] = 1.16
        paramsGaAs["c11"] = 1.223
        paramsGaAs["c12"] = 0.571
        paramsGaAs["c44"] = 0.600
        paramsGaAs["vlong"] = 4730
        paramsGaAs["massdens"] = 5317
        paramsGaAs["molV"] = 0.04517
        paramsGaAs["lattconst"] = 5.653
        super(GaAs,self).__init__(name,paramsGaAs)
        
    def copy(self):
        return GaAs(self.name)

class AlAs(Material):
    '''AlAs, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe]. Elastic constants and hydrostatic potentials from 
    [vandeWalle1989].
    '''
    def __init__(self,name = None):
        if name is None:
            name = "AlAs"
        paramsAlAs = Material.params_dict.copy()
        
        paramsAlAs["meff"] = 0.15
        paramsAlAs["Ec"] = 0.99
        paramsAlAs["Eg"] = 3.099
        paramsAlAs["EX"] = 2.24
        paramsAlAs["ELO"] = 0.050
        paramsAlAs["Ep"] = 22
        paramsAlAs["eps0"] = 10.06
        paramsAlAs["epsinf"] = 8.05
        paramsAlAs["lattconst"] = 5.6611
        paramsAlAs["ac"] = -7.17
        paramsAlAs["av"] = 1.16
        paramsAlAs["c11"] = 1.250
        paramsAlAs["c12"] = 0.534
        paramsAlAs["c44"] = 0.542
        paramsAlAs["molV"] = 0.04536
        paramsAlAs["massdens"] = 3760
        paramsAlAs["vlong"] = 5650
        super(AlAs,self).__init__(name,paramsAlAs)
        
    def copy(self):
        return AlAs(self.name)

class InAs(Material):
    '''InAs, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe]. Elastic constants and hydrostatic potentials from 
    [vandeWalle1989].
    '''    
    def __init__(self,name = None):
        
        if name is None:
            name = "InAs"
            params = Material.params_dict.copy()
            params["meff"] = 0.026
            params["Ec"] = -0.892
            params["Eg"] = 0.417
            params["EX"] = 1.433
            params["ELO"] = 0.0298 # Ioffe
            params["Ep"] = 21.5
            params["eps0"] = 15.1 # Ioffe
            params["epsinf"] = 12.3
            params["ac"] = -5.08
            params["av"] = 1.00
            params["c11"] = 0.833
            params["c12"] = 0.453
            params["c44"] = 0.396
            params["vlong"] = 3830
            params["lattconst"] = 6.0583
            params["massdens"] = 5680
            params["molV"] = 0.0556
            
            super(InAs,self).__init__(name,params)
            
    def copy(self):
        return InAs(self.name)
    
class InP(Material):
    '''InP, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe]. Elastic constants and hydrostatic potentials from 
    [vandeWalle1989].
    '''    
    def __init__(self,name = None):
        
        if name is None:
            name = "InP"
            params = Material.params_dict.copy()
            params["meff"] = 0.0795
            params["Ec"] = -0.2854 # Vurgaftman VBO GaAs/InP = 0.19
            params["Eg"] = 1.4236
            params["ELO"] = 0.043 # Ioffe
            params["Ep"] = 20.7
            params["eps0"] = 12.5 # Ioffe
            params["epsinf"] = 9.61
            params["ac"] = -5.04
            params["av"] = 1.27
            params["c11"] = 1.022
            params["c12"] = 0.576
            params["c44"] = 0.460
            params["vlong"] = 4580
            params["lattconst"] = 5.8687
            params["massdens"] = 4810
            params["molV"] = 0.05056
            super(InP,self).__init__(name,params)
            
    def copy(self):
        return InP(self.name)
    
class GaSb(Material):
    '''GaSb, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe]. Elastic constants and hydrostatic potentials from 
    [vandeWalle1989].
    '''    
    def __init__(self,name = None):
        
        if name is None:
            name = "GaSb"
            params = Material.params_dict.copy()
            params["meff"] = 0.039
            params["Ec"] = 0.063
            params["Eg"] = 0.812
            params["ELO"] = 0.0297 # Ioffe
            params["Ep"] = 27.0
            params["eps0"] = 15.7 # Ioffe
            params["epsinf"] = 14.4
            params["ac"] = -7.5
            params["av"] = 0.79
            params["c11"] = 0.908
            params["c12"] = 0.413
            params["c44"] = 0.445
            params["vlong"] = 3970
            params["lattconst"] = 6.0959
            params["massdens"] = 5610
            params["molV"] = 0.05663
            super(GaSb,self).__init__(name,params)
            
    def copy(self):
        return GaSb(self.name)
    
class AlSb(Material):
    '''GaSb, band parameters from [Vurgaftman2001]. Other parameters from
    [Adachi1999]. Elastic constants and hydrostatic potentials from 
    [vandeWalle1989].
    '''    
    def __init__(self,name = None):
        
        if name is None:
            name = "AlSb"
            params = Material.params_dict.copy()
            params["meff"] = 0.14
            params["Ec"] = 1.257
            params["Eg"] = 2.386 #direct gap (indirect gap =1.6 eV)
            params["ELO"] = 0.0298
            params["Ep"] = 18.7
            params["eps0"] = 12.04
            params["epsinf"] = 10.24
            params["ac"] = -4.5
            params["av"] = 1.38
            params["c11"] = 0.877
            params["c12"] = 0.434
            params["c44"] = 0.408
            params["vlong"] = 3830
            params["lattconst"] = 6.1355
            params["massdens"] = 5680
            params["molV"] = 0.05774
            super(AlSb,self).__init__(name,params)
            
    def copy(self):
        return AlSb(self.name)
    
class InSb(Material):
    '''GaSb, band parameters from [Vurgaftman2001]. Other parameters from
    [Ioffe]. Elastic constants and hydrostatic potentials from 
    [vandeWalle1989].
    '''    
    def __init__(self,name = None):
        
        if name is None:
            name = "InSb"
            params = Material.params_dict.copy()
            params["meff"] = 0.0135
            params["Ec"] = -0.484
            params["Eg"] = 0.235
            params["ELO"] = 0.025 # Ioffe
            params["Ep"] = 23.3
            params["eps0"] = 16.8 # Ioffe
            params["epsinf"] = 15.7
            params["ac"] = -6.94
            params["av"] = 0.36
            params["c11"] = 0.659
            params["c12"] = 0.356
            params["c44"] = 0.300
            params["vlong"] = 3400
            params["lattconst"] = 6.4794
            params["massdens"] = 5770
            params["molV"] = 0.06800
            super(InSb,self).__init__(name,params)
            
    def copy(self):
        return InSb(self.name)

class ZnO(Material):
    '''ZnO, parameters from [Janotti2007] and [Bateman1962].
    Note that at this time, the state of knowledge about this
    material is very limited.'''
    def __init__(self,name = None):
        if name is None:
            name = "ZnO"
        params = Material.params_dict.copy()
        params["meff"] = 0.22
        params["Ec"] = -0.44 # ±0.23 [ZhangAPL2008]
        params["Eg"] = 3.4 # [Janotti2007]
        params["ELO"] = 0.072
        params["Ep"] = 21.5
        params["eps0"] = 8.49
        params["epsinf"] = 3.72
        params["ac"] = -2.3 # [Janotti2007]
        params["av"] = -0.6 # [Janotti2007]
        params["c11"] = 0.419 # JemmyCinitA 2014
        params["c12"] = 0.0384 # JemmyCinitA 2014
        params["c44"] = -0.0587 # JemmyCinitA 2014
        params["vlong"] = 6090 # Bateman JAP 33, 1962
        params["massdens"] = 5606
        params["molV"] = 0.0145
        params["lattconst"] = 5.2
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
        params = Material.params_dict.copy()
        params["meff"] = 0.37
        params["Ec"] = 1.26 # Janotti, van der Walle PRB 2007 (-0.44 eV)
        params["Eg"] = 6.3 # Janotti, van der Walle PRB 2007
        params["ELO"] = 0.089
        params["Ep"] = 21.5
        params["eps0"] = 9.6
        params["epsinf"] = 2.98
        params["ac"] = -4.3 # [Janotti2007]
        params["av"] =  2.0 # [Janotti2007]
        params["c11"] = 0.2979 # Strauch, "Semiconductors" 2017
        params["c12"] = 0.0958 # Strauch, "Semiconductors" 2017
        params["c44"] = 0.1542 # Strauch, "Semiconductors" 2017
        params["lattconst"] = 4.21
        super(MgO,self).__init__(name,params)
        
    def copy(self):
        return MgO(self.name)
    
class MgOzoterac(Material):
    '''MgO, parameters used in the Zoterac (ERC) project.'''
    
    def __init__(self,name = 'MgOzoterac'):
        params = Material.params_dict.copy()
        params["meff"] = 0.22
        params["Ec"] = 0.86 # 1.3-0.44 Deliverable D1.3 1/9/2017
        params["Eg"] = 5.3 # Deliverable D1.3 1/9/2017
        params["ELO"] = 0.089
        params["Ep"] = 21.5
        params["eps0"] = 9.6
        params["epsinf"] = 2.98
        params["ac"] = -4.3
        params["lattconst"] = 4.21
        super(MgOzoterac,self).__init__(name,params)

    def copy(self):
        return MgOzoterac(self.name)
        
# Ternaries:

class AlGaAs(Material):
    '''Al_xGa_1-xAs. Bowing parameters from [Vurgaftman2001]. 
    CBO from Yi et al, PRB 81, 235325 (2010)'''

    def __init__(self, name = None, x=0., CBO_Yi = True):
        self.CBO_Yi = CBO_Yi
        if name is None:
            name = "Al_" + str(x) + "GaAs"
        mat1 = AlAs()
        mat2 = GaAs()
        C = Material.params_dict.copy()
        C["Eg"] = -0.127 + 1.310*x
        C["EX"] = 0.055
        C["Ec"] = -0.127 + 1.310*x
        super(AlGaAs,self).__init__(name, Material.params_dict.copy(), 
                                    mat1, mat2, C, x, GaAs())
            
    def updateAlloy(self,x,reset_strain = False):
        self.C["Eg"] = -0.127 + 1.310*x
        self.C["Ec"] = -0.127 + 1.310*x
        self.C["EX"] = 0.055
        super(AlGaAs,self).updateAlloy(x,reset_strain = reset_strain)
        if self.CBO_Yi:
            if x < 0.42:
                self.params["Ec"] = 0.831*x
            else:
                #self.params["Ec"] = 0.332 + 0.054*x #indirect gap
                self.params["Ec"] = -0.115 + 1.105*x  #direct gap
            
    def copy(self):
        return AlGaAs(self.name,self.x)
    

class InGaAs_on_GaAs(Material):
    '''
    In_xGa_1-xAs on GaAs. Bowing parameters from [ArentJAP1989].
    Tested experimentally up to x = 0.28.
    '''
    def __init__(self,name = None, x = 0.):
        if name is None:
            name = "In_" + str(x) + "GaAs/GaAs"
        mat1 = InAs()
        mat2 = GaAs()
        mat1.params["Eg"] = 0.435
        mat1.params["Ec"] = -0.8672
        A = Material.params_dict.copy()
        
        A["Eg"] = -0.496
        A["meff"] = 0.0091 # same as InP
        A["Ep"] = -1.48 # same as InP
        A["ac"] = 2.61 # same as InP
        A["eps0"] = -0.67 # same as InP
        A["ELO"] = 0.002 # same as InP
        A["Ec"] = 0.397
        super(InGaAs_on_GaAs,self).__init__(name,Material.params_dict.copy(),
                                            mat1, mat2, A, x, GaAs())
        
    def copy(self):
        return InGaAs_on_GaAs(self.name,self.x)
    
    
    
class InGaAs(Material):
    '''In_xGa_1-xAs on InP. Bowing parameters from [Vurgaftman2001] and [Ioffe].'''
    def __init__(self,name = None, x = None):
        if x is None:
            x = 0.53
        if name is None:
            name = "In_" + str(x) + "GaAs"
        mat1 = InAs()
        mat2 = GaAs()
        A = Material.params_dict.copy()
        
        A["Eg"] = 0.477
        A["meff"] = 0.0091
        A["Ep"] = -1.48
        A["ac"] = 2.61
        A["eps0"] = -0.67 # Ioffe
        A["ELO"] = 0.002 # fit to Ioffe 34 meV lattice matched
        A["Ec"] = 0.097 # (= Cvbo + Cgap, Vurgaftman)
        super(InGaAs,self).__init__(name,Material.params_dict.copy(),
                                    mat1, mat2, A, x, InP())
        
    def copy(self):
        return InGaAs(self.name,self.x)
    

        
class AlInAs(Material):
    '''Al_xIn_1-xAs on InP. Bowing parameters from [Vurgaftman2001].'''
    def __init__(self,name = None, x = None):
        if x is None:
            x = 0.48
        if name is None:
            name = "Al_" + str(x) + "InAs"
        
        mat1 = AlAs()
        mat2 = InAs()
        A = Material.params_dict.copy()
        A["Eg"] = 0.70
        A["meff"] = 0.049
        A["Ep"] = -4.81
        A["ac"] = -1.4
        A["Ec"] = -0.0383 # (= Cvbo + Cgap + mod, Vurgaftman, mod to match LM DeltaEc)
     
            
        super(AlInAs,self).__init__(name,Material.params_dict.copy(),
                                    mat1, mat2, A, x, InP())
        
    def copy(self):
        return AlInAs(self.name,self.x)
        
class GaInSb(Material):
    '''Ga_xIn_(1-x)Sb/GaSb. Bowing parameters from [Vurgaftman2001].'''
    
    def __init__(self, name = None, x = None):
        if x is None:
            x = 0.0
        if name is None:
            name = "Ga_" + str(x) + "InSb"
        
        mat1 = GaSb()
        mat2 = InSb()
        A = Material.params_dict.copy()
        A["Eg"] = 0.415
        A["EX"] = 0.33
        A["meff"] = 0.0092
        A["Ec"] = 0.415 # (= Cvbo + Cgap, Vurgaftman)
     
            
        super(GaInSb,self).__init__(name,Material.params_dict.copy(),
                                    mat1, mat2, A, x, GaSb())
        
    def copy(self):
        return GaInSb(self.name,self.x)

class AlInSb(Material):
    '''Al_xIn_(1-x)Sb/GaSb. Bowing parameters from [Vurgaftman2001].'''
    
    def __init__(self, name = None, x = None):
        if x is None:
            x = 0.0
        if name is None:
            name = "Al_" + str(x) + "InSb"
        
        mat1 = AlSb()
        mat2 = InSb()
        A = Material.params_dict.copy()
        A["Eg"] = 0.43
        A["Ec"] = 0.43 # (= Cvbo + Cgap, Vurgaftman)
            
        super(AlInSb,self).__init__(name,Material.params_dict.copy(),
                                    mat1, mat2, A, x, GaSb())
        
    def copy(self):
        return AlInSb(self.name,self.x)

class AlGaSb(Material):
    '''Al_xGa_(1-x)Sb/GaSb. Bowing parameters from [Vurgaftman2001].'''
    
    def __init__(self, name = None, x = None):
        if x is None:
            x = 0.0
        if name is None:
            name = "Al_" + str(x) + "InSb"
        
        mat1 = AlSb()
        mat2 = GaSb()
        A = Material.params_dict.copy()
        A["Eg"] = -0.044 + 1.22*x
        A["Ec"] = -0.044 + 1.22*x # (= Cvbo + Cgap, Vurgaftman)
            
        super(AlGaSb,self).__init__(name,Material.params_dict.copy(),
                                    mat1, mat2, A, x, GaSb())
    def updateAlloy(self,x,reset_strain = False):
        
        self.C["Eg"] = -0.044 + 1.22*x 
        self.C["Ec"] = -0.044 + 1.22*x  
        
        super(AlGaSb,self).updateAlloy(x,reset_strain = reset_strain)
               
     
    def copy(self):
        return AlGaSb(self.name,self.x)

class GaAsSb(Material):
    '''Ga_xAs_(1-x)Sb/GaSb. Bowing parameters from [Vurgaftman2001].'''
    
    def __init__(self, name = None, x = None):
        if x is None:
            x = 0.0
        if name is None:
            name = "Ga_" + str(x) + "AsSb"
        
        mat1 = GaAs()
        mat2 = GaSb()
        A = Material.params_dict.copy()
        A["Eg"] = 1.43
        A["Ec"] = 0.37 # (= Cvbo + Cgap, Vurgaftman)
            
        super(GaAsSb,self).__init__(name,Material.params_dict.copy(),
                                    mat1, mat2, A, x, GaSb())
        
    def copy(self):
        return GaAsSb(self.name,self.x)
    
class InAsSb(Material):
    '''In_xAs_(1-x)Sb/GaSb. Bowing parameters from [Vurgaftman2001].'''
    
    def __init__(self, name = None, x = None):
        if x is None:
            x = 0.0
        if name is None:
            name = "InAs_" + str(x) + "Sb"
        
        mat1 = InAs()
        mat2 = InSb()
        A = Material.params_dict.copy()
        A["Eg"] = 0.67
        A["meff"] = 0.035
        A["Ec"] = 0.67 # (= Cvbo + Cgap, Vurgaftman)
            
        super(InAsSb,self).__init__(name,Material.params_dict.copy(),
                                    mat1, mat2, A, x, GaSb())
        
    def copy(self):
        return InAsSb(self.name,self.x)
    
class AlAsSb(Material):
    '''AlAs_(1-x)Sb_x/GaSb. Bowing parameters from [Vurgaftman2001].'''
    
    def __init__(self, name = None, x = None):
        if x is None:
            x = 0.0
        if name is None:
            name = "AlAs_(1-" + str(x) + ")Sb_x"
        
        mat1 = AlSb()
        mat2 = AlAs()
        A = Material.params_dict.copy()
        A["Eg"] = 0.8
        A["Ec"] = -0.91 # (= Cvbo + Cgap, Vurgaftman)
            
        super(AlAsSb,self).__init__(name,Material.params_dict.copy(),
                                    mat1, mat2, A, x, GaSb())
        
    def copy(self):
        return AlAsSb(self.name,self.x)
    
class ZnMgO(Material):
    '''Zn_xMg_1-xO. Only linear mixing.'''
    def __init__(self,name = None,x = 0.):
        if name is None:
            name = "Zn_"+str(x)+"Mg_"+str(1-x)+"O"
        mat1 = ZnO()
        mat2 = MgO()
        # Bowing paremeters:Unknown
        A = Material.params_dict.copy()
        super(ZnMgO,self).__init__(name,Material.params_dict.copy(),mat1,mat2,A,x,ZnO())
        
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
                print("General case is not implemented, use y = 0.48 or None!")
                
            if z is None:
                z = 0.53
            elif z != 0.53:
                print("General case is not implemented, use z = 0.53 or None!")
            
            mat1 = AlInAs()
            mat2 = InGaAs()
            
            A = Material.params_dict.copy()
            A["Eg"] = 0.22
            A["meff"] = -0.016
            super(AlInGaAs,self).__init__(name,Material.params_dict.copy(),mat1,mat2,A,x)
            
        def copy(self):
            return AlInGaAs(self.x, self.name)
        
        def updateAlloy(self, x, reset_strain = False):
            Material.updateAlloy(self, x, reset_strain = reset_strain)
            self.params["Valloy"] = self.mat1.params["Valloy"]*x + \
                self.mat2.params["Valloy"]*(1-x)
            self.params["Ec"] = 0.73*(0.712*x - 0.22*x*(1-x)) + self.mat2.params["Ec"]
        