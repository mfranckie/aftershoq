'''
Created on 26 Jan 2018

@author: martin
'''
import random
import hilbert
import numpy as np

class MaterialPar():
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
        self.material = material
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
        if orig is None:
            self.Nl = 0;
            self.layers = []
            self.dopings = []
        else:
            self.Nl = 0
            self.layers = []
            for l in orig.layers:
                self.addLayer(Layer(l.width,l.material,l.eta,l.lam))
            self.dopings = orig.dopings[:]
        
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

class Material:
    def __init__(self,name,params,mat1=None,mat2=None,C=None,x=None):
        self.name = name
        self.params = params
        if mat1 is not None:
            self.mat1 = mat1
        if mat2 is not None:
            self.mat2 = mat2
        if C is not None:
            self.C = C
        if x is not None:
            self.x = x
            self.updateAlloy(x)
            
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
        
    
class Sgenerator():
    def __init__(self,origstruct,dw,dx,ddop=None):
        self.dx = dx
        self.dw = dw
        self.orig = origstruct
        self.structures = []
        if ddop is None:
            ddop = [0,0,0]
        self.ddop = ddop
        for l in origstruct.layers:
            if l.width-dw[origstruct.layers.index(l)] <= 0:
                dw[origstruct.layers.index(l)] = l.width-0.01
    
    def genRanStructs(self,N):
        """
        Generate N structures with random variations in layer widths
        """
        for _ in range(0,N):
            news = Structure(self.orig)
            
            for il in range(0,len(self.dw)):
                w0 = self.orig.layers[il].width
                width = w0 + 2*(random.random()-0.5)*self.dw[il]
                news.layers[il].width = width
                
            news.dopings = []
            for dl in self.orig.dopings:
                
                # find out offset from layer start for original structure
                zc = (dl[0]+dl[1])/2
                layerindex = self.orig.layerIndex(zc)
                dzc = zc - self.orig.layerPos(layerindex)
                               
                # define new offset from layer start
                dzc = dzc + 2*(random.random()-0.5)*self.ddop[0]
                dopw = (dl[1]-dl[0]) + 2*(random.random()-0.5)*self.ddop[1]
                dopdens = dl[2] + 2*(random.random()-0.5)*self.ddop[2]
                
                # add new doping layer
                zi = dzc - dopw/2
                zf = dzc + dopw/2
                news.addDoping(zi, zf, dopdens, layerindex)
                
            self.structures.append(news)
            del news
           
    def genRanHilbertStructs(self,N,p):
        """
        Generate N structures on the hilbert curve defined by p
        Number of dimensions is determined by # parameters to 
        change
        """
        
        # Number of dimensions in Euclidian space
        ND = 0
        
        # Defince list of changeable parameters
        params = []
        dparams = []
        windex = []
        xindex = []
        dopindex = []
        for i in range(0,len(self.orig.layers)):
            if self.dw[i]>0:
                params.append(self.orig.layers[i].width)
                dparams.append(self.dw[i])
                windex.append(i)
                ND += 1
        # TODO: extend to multiple doping layers
        for i in range(0,len(self.orig.dopings[0])):
            if self.ddop[i] > 0:
                params.append(self.orig.dopoings[0][i])
                dparams.append(self.ddop[i])
                dopindex.append(i)
                ND +=1
        for i in range(0,len(self.dx)):
            if self.dx[i] > 0:
                params.append(self.orig.layers[i].material.x)
                dparams.append(self.dx[i])
                xindex.append(i)
                ND += 1
  
        
        # define hilbert curve:
        hilbert_curve = hilbert.HilbertCurve(p, ND)
        
        # imax = length of hilbert curve (# nodes)
        imax = 2**(ND*p)-1
        
        # pmax = side of hyper cube
        pmax = 2**p-1
        
        print 'Dim = {0}, p = {1}, imax = {2}, pmax = {3}'.format(ND,p,imax,pmax)
        
        #define distance along hilbert curve:
        dhil = []
        [dhil.append(i) for i in range (0,imax+1) ]
        
        # define mapping from distance to coordinates
        hmap = []
        [hmap.append(hilbert_curve.coordinates_from_distance(i)) for i in dhil]
        
        # scaling each dimension to p_iE[0:pmax]:
        
        arr_par = np.array(params)
        arr_dpar = np.array(dparams)
        parmin = arr_par - arr_dpar
        parmax = arr_par + arr_dpar
        
        # contains all parameters, scaled to the range [0:pmax]
        par_scaled = (arr_par - parmin)/(parmax-parmin)*pmax
        
        hmap_array = np.array(hmap)
        
        coordinates = []
        w1 = []
        w2 = []
        for _ in range(0,N):
            # pick a random point along the hilbert curve:
            d = random.random()*imax
            
            # interpolate to get the coordingates:
            intd = int(d)
            dd = float(dhil[intd])-d
            
            if intd == imax:
                dp = hmap_array[intd-1][:]-hmap_array[intd][:]
            else:
                dp = hmap_array[intd][:]-hmap_array[intd+1][:]
            
            # x_ipl now contains the scaled coordingates (widths, x:s, and dopings)
            p_ipl = hmap_array[intd][:] + dp*dd
            
            # un-scale it and create structure
            p_unscaled = p_ipl*(parmax-parmin)/pmax + parmin
            
            # generate strucure:
            news = Structure(self.orig)
            pindex = 0
            for i in range(0,len(windex)):
                news.layers[windex[i]].width = p_unscaled[pindex]
                pindex+=1
            for i in range(0,len(dopindex)):
                news.dopings[dopindex[i]] = p_unscaled[pindex]
                pindex +=1
            for i in range(0,len(xindex)):
                news.layers[xindex[i]].material.updateAlloy(p_unscaled[pindex])
                pindex +=1
            self.structures.append(news)
            #del news
            coordinates.append(p_unscaled)
            w1.append(p_unscaled[0])
            w2.append(p_unscaled[1])
       
        #pl.plot(w1,w2,'*')
        x1 = []
        y1 = []
        [x1.append(float(hmap[i][0])) for i in range(0,imax+1)]
        [y1.append(float(hmap[i][1])) for i in range(0,imax+1)]
        x1 = np.array(x1)/pmax*(parmax[0]-parmin[0]) + parmin[0]
        y1 = np.array(y1)/pmax*(parmax[1]-parmin[1]) + parmin[1]
        #pl.plot(x1,y1,'-+',hold = 'true')
        #pl.show()
        return coordinates
