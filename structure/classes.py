'''
Created on 26 Jan 2018

@author: martin

Module containing core classes for materials, layers, and structures.
'''

import structure.matpar as mp

class Layer:
    ''' Defines a layer in the heterostructure. Each layer has a material
    or alloy, interface roughness mean height eta and correlation length
    lam. These are defined as the roughness parameters for the interface
    between this one and the following layer.
    
    ''' 
    
    def __init__(self,width,material,eta,lam):
        '''Constructor. Parameters:
        width:  Width of the layer in nm
        material:  Material object
        eta:  Interface roughness mean height
        lam:  Interface roughness correlation length
        '''
        
        self.width = width
        self.material = material.copy()
        self.lam = lam
        self.eta = eta
        
    def __str__(self):
        return str([self.width,self.material,self.eta,self.lam])
        
    def __repr__(self):
        return str(self)


class Structure:
    '''Defines an entire heterostructure, consisting of several
    Layers and doping regions.
    '''
    
    sid = 0
    
    def __init__(self, orig=None):
        ''' Constructur. Optionally copies from Structure object orig.
        Each structure receives an id (sid) and a dirname, which defaults to 
        sid.
        '''
        
        self.sid = Structure.sid
        Structure.sid +=1
        
        self.dirname = str(self.sid)
        self.length = 0
        self.dopings = []
        if orig is None:
            self.Nl = 0;
            self.layers = []
        else:
            self.Nl = 0
            self.layers = []
            for l in orig.layers:
                self.addLayer(Layer(l.width,l.material,l.eta,l.lam))
            for dl in orig.dopings:
                self.addDoping(dl[0], dl[1], dl[2])
        
    def setIFR(self,eta,lam):
        ''' Manually set interface roughness parameters.'''
        
        self.eta = eta
        self.lam = lam
    
    def addLayer(self,layer):
        '''Add a Layer to this structure.'''
        self.layers.append(layer)
        self.Nl+=1
        self.length += layer.width
    
    def addLayerIFR(self,width,material,eta,lam):
        '''Add a new Layer, created from 
        width:   Width in nm
        material:  Material object for new layer
        eta:  Interface roughness mean height
        lam:  Interface roughness correlation length
        
        '''
        self.addLayer(Layer(width,material,eta,lam))
        
    def addLayerMW(self,width,material):
        '''Add a new layer based on width (nm) and Material.
        Uses pre-defined interface roughness parameters (set via setIFR )
        '''
        self.addLayer(Layer(width,material,self.eta,self.lam))
        
    def addDoping(self,zi,zf,density,layerindex = None):
        '''Add a doping layer.
        zi:  starting position of layer
        zf:  end position of layer
        density:  volume doping density (cm^-3)
        layerindex: optional, set (zi, zf) relative to Layer start
        '''
        
        if layerindex is not None:
            lp = self.layerPos(layerindex)
            zi += lp
            zf += lp
        self.dopings.append([zi,zf,density])
        
    def layerPos(self,index):
        '''Returns the position of Layer with index "index"
        '''
        pos = 0
        for l in self.layers[0:index]:
            pos += l.width
        return pos
    
    def layerIndex(self,pos):
        '''Returns the index of the Layer covering the position pos.'''
        
        z = 0
        for li in self.layers:
            z += li.width
            if pos<z:
                return self.layers.index(li)
    
    def __str__(self):
        return str(self.layers)
    
    def layername(self):
        '''Generates and returns a name for this structure, based on 
        layer widths.
        '''
        name = ""
        for l in self.layers:
            name = name + str(l.width) + "_"
        return name
        
    layers=[]
    dopings=[]

class Material(object):
    '''Defines a material or alloy between two materials.'''
    
    def __init__(self,name,params,mat1=None,mat2=None,C=None,x=None):
        '''Constructor. Parameters:
        name:  The name of this material
        params:  Material parameters. Set via module matpar.
        mat1:   Optional, Material one in new alloy
        mat2:   Optional, Material two in new alloy
        C:      Optional, Bowing parameters for new alloy. All default to 0.
        x:      Optional, relative compositon of mat1 to mat2
        '''
        self.name = name
        self.params = params
        self.mat1 = mat1
        self.mat2 = mat2
        self.C = C
        self.x = x
        if x is not None:
            self.updateAlloy(x)
            
    def copy(self):
        ''' Returns a deep copy of this material.'''
        return Material(self.name, self.params, self.mat1, self.mat2, self.C, self.x)
            
    def updateAlloy(self,x):
        '''Updates the alloy composition of this alloy material with x
        as the new composition.
        '''
        
        if(self.x is not None):
            self.x = x
            self.params = self.alloy(self.mat1,self.mat2,self.C,self.x)
        
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return str(self)
    
    @staticmethod
    def alloy(mat1,mat2,C,x):
        '''Creates and returns a new set of material parameters, which are
        defined by the parameters of mat1, mat2, the bowing parameters
        in C, and the composition x. C is a list of bowing parameters matching
        the parameters in Material.params and MaterialPar:
        newparams = x*param1 + (1-x)*param2 - C*x*(1-x)
        x is the alloy fraction of mat1 to mat2.
        '''
        
        params3=[]
        for i in range(0,mp.Nparam):
            params3.append(mat1.params[i]*x + mat2.params[i]*(1-x) - C[i]*x*(1-x))
        # Alloy scattering:
        params3[mp.Valloy] = x*(1-x)*(mat1.params[mp.Ec]-mat2.params[mp.Ec])
        return params3
