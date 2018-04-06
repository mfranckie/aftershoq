'''
Created on 26 Jan 2018

@author: martin
'''

class MaterialPar:
    meff = 0
    Ec = 1
    Eg = 2
    Ep = 3
    Valloy = 4
    ELO = 5
    eps0 = 6
    epsinf = 7
    Vdef = 8
    vlong = 9
    massdens = 10
    molV = 11
    lattconst = 12
    Nparam = 13
    
    valdict = {"meff":meff, "CBO": Ec, "Eg": Eg, "Ep":Ep, "Alloy pot.": Valloy,
               "ELO": ELO, "eps(0)":eps0, "eps(inf)":epsinf,
               "deform. pot.":Vdef,"long. sound vel.":vlong,"mass dens":massdens,
               "mol volume":molV,"lattice constant": lattconst}
    
    paramList = []
    
    def __init__(self):
        self.initList(self.paramList)
    
    def setParam(self,param,value):
        self.paramList[param] = value
    
    @classmethod
    def initList(cls,C):
        for _ in range(0,cls.Nparam):
            C.append(0)

class Layer:
    def __init__(self,width,material,eta,lam):
        self.width = width
        self.material = material.copy()
        self.lam = lam
        self.eta = eta
        
    def __str__(self):
        return str([self.width,self.material,self.eta,self.lam])
        
    def __repr__(self):
        return str(self)

class Structure:
    
    sid = -1
    
    def __init__(self, orig=None):
        self.sid = Structure.sid
        Structure.sid +=1
        # default dirname
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
 
    def copy(self,structure):
        self.layers = []
        for l in range(0,structure.Nl):
            self.addLayer(structure.layers[l])
        
    def setIFR(self,eta,lam):
        self.eta = eta
        self.lam = lam
    
    def addLayer(self,layer):
        self.layers.append(layer)
        self.Nl+=1
        self.length += layer.width
    
    def addLayerIFR(self,width,material,eta,lam):
        self.addLayer(Layer(width,material,eta,lam))
        
    def addLayerMW(self,width,material):
        self.addLayer(Layer(width,material,self.eta,self.lam))
        
    def addDoping(self,zi,zf,density,layerindex = None):
        if layerindex is not None:
            lp = self.layerPos(layerindex)
            zi += lp
            zf += lp
        self.dopings.append([zi,zf,density])
        
    def layerPos(self,index):
        pos = 0
        for l in self.layers[0:index]:
            pos += l.width
        return pos
    
    def layerIndex(self,pos):
            z = 0
            for li in self.layers:
                z += li.width
                if pos<z:
                    return self.layers.index(li)
    
    def __str__(self):
        return str(self.layers)
    
    def layername(self):
        name = ""
        for l in self.layers:
            name = name + str(l.width) + "_"
        return name
        
    layers=[]
    dopings=[]

class Material(object):
    def __init__(self,name,params,mat1=None,mat2=None,C=None,x=None):
        self.name = name
        self.params = params
        self.mat1 = mat1
        self.mat2 = mat2
        self.C = C
        self.x = x
        if x is not None:
            self.updateAlloy(x)
            
    def copy(self):
        return Material(self.name, self.params, self.mat1, self.mat2, self.C, self.x)
            
    def updateAlloy(self,x):
        if(self.x is not None):
            self.x = x
            self.params = self.alloy(self.mat1,self.mat2,self.C,self.x)
        
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return str(self)
    
    @staticmethod
    def alloy(mat1,mat2,C,x):
        params3=[]
        for i in range(0,MaterialPar.Nparam):
            params3.append(mat1.params[i]*x + mat2.params[i]*(1-x) - C[i]*x*(1-x))
        # Alloy scattering:
        params3[MaterialPar.Valloy] = x*(1-x)*(mat1.params[MaterialPar.Ec]-mat2.params[MaterialPar.Ec])
        return params3
