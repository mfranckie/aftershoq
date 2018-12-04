'''
Created on 26 Jan 2018

@author: martin

Module containing core classes for materials, layers, and structures.
'''

from . import matpar as mp
import copy
import numpy as np

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
        self.material = copy.deepcopy(material)
        
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
    
    def __init__(self, orig=None, name = None):
        ''' Constructur. Optionally copies from Structure object orig.
        Each structure receives an id (sid) and a dirname, which defaults to 
        sid.
        '''
        
        self.sid = Structure.sid
        Structure.sid +=1
        
        if name is None:
            self.name = type(self).__name__
        else:
            self.name = name
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
        
        for l in self.layers:
            l.eta = self.eta
            l.lam = self.lam
    
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
        
    def addLayerWM(self,width,material):
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
            
    def layerDoping3D(self, index):
        ''' Returns the sheet doping density in the layer with layer index
        "Index" (in units of cm^-2*nm^-1).'''
        doping = 0
        l0 = self.layerPos(index)
        l1 = self.layerPos(index) + self.layers[index].width
        for dop in self.dopings:
            if dop[1] >= l0 and dop[0] <= l1:
                if dop[0] > l0:
                    if dop[1] > l1:
                        ol = l1 - dop[0]
                    else:
                        ol = dop[1]-dop[0]
                else:
                    if dop[1] > l1:
                        ol = l1 - l0
                    else:
                        ol = dop[1]-l0
                    
                    
                doping += dop[2]*ol/(dop[1]-dop[0])
                
        return doping
    
    def convert_to_ML(self):
        """
        Converts this structure to integer monolayers. Assumes doping
        spans entire layers and keeps sheet density constant.
        Lattice constants are given in Å, layer widths in nm.
        """
        
        dops = []
        [dops.append(self.layerDoping3D(i)) for i in range(len(self.layers))]
        print("dops = ", dops)
        
        self.dopings = []
        
        for i in range( len(self.layers) ):
            layer = self.layers[i]
            NML = np.round(layer.width/layer.material.params[mp.lattconst]*10.)
            layer.width = NML*layer.material.params[mp.lattconst]/10.
            if(dops[i]>0):
                self.addDoping(0, layer.width, dops[i], i)
            
    
    def __str__(self):
        s = "[width, Material, eta, lambda] (id="+str(self.sid) + ")\n" 
        for l in self.layers:
            s += str(l) + "\n"
        return s
    
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
    
    def getSheetDop(self):
        '''Returns the sheet doping density per period in cm^-2.'''
        sd = 0.
        for l in self.dopings:
            sd += (l[1]-l[0])*l[2]
        return sd*1e-7
    
    def prettyPrint(self):
        s = self.name + ":\n"
        s += "width\tMaterial\tVol. doping\n(nm)\t\t\t(cm-3)\n"
        for i in range(0, len(self.layers)):
            l = self.layers[i]
            s+= str(l.width)+"\t"+str(l.material)+"\t"+str(self.layerDoping3D(i))+"\n"
        return s

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
        if (x is not None):
            if x > 1:
                print("ERROR: x > 1 in material creation! Stopping.")
                exit(1)
            self.updateAlloy(x)
            
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
