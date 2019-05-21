'''
Created on 26 Jan 2018

@author: martin

Module containing core classes for materials, layers, and structures.
'''

import copy
import numpy as np
import os.path as Path
from lxml import etree
import datetime
import scipy.optimize as opt

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

    postfix = ".aqst"

    sid = 0

    def __init__(self, orig=None, name = None, T = 0):
        """Constructor. Optionally copies from Structure object orig.
        Each structure receives an id (sid) and a dirname, which defaults to
        sid.

        Parameters
        ----------
        orig : Structure
            (Optional) Structure to copy from.
        name : String
            (Optional) Name of the structure.
        T : float
            (Optional) Lattice temperature that should be used for every
            material used in this structure.
        """

        self.TL = T
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

    def calcStrain(self):
        """
        Calculates return the total strain in one period h, and the averaged
        critical strain per period Sc. The structure is estimated to relax when
        abs(h) > abs(Sc) (or after N periods where N * abs(h) > abs(Sc) ).

        Returns: h, Sc
        """

        # All layers must have the same substrate
        asub = self.layers[0].material.substrate.params["lattconst"]

        h = 0. # total strain per period, to be summed
        Sc = 0. # critical strain averaged over two periods
        for l in self.layers*2:
            l.material.calcStrain()
            h += (l.material.params["lattconst"] - asub)/asub*l.width
            Sc += (l.material.params["lattconst"] - asub)/asub*l.material.hcrit()

        h /= 2.
        # averaged critical strain per period (absolute value)
        Sc /= len(self.layers)*2
        return h, Sc


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
            NML = np.round(layer.width/layer.material.params["lattconst"]*10.*2.)
            layer.width = NML*layer.material.params["lattconst"]/10./2.
            if(dops[i]>0):
                self.addDoping(0, layer.width, dops[i], i)

    def compensate_strain(self, x=None, y=None, Ec = None):
        """Compensate the stain within one period by adjusting the alloy
        compositions. Attempts to keep the band offset the same.

        Parameters
        ----------
        x : float
            (Optional) Initial alloy composition of layers[0].material.
        y : float
            (Optional) Initial alloy composition of layers[1].material.
        Ec : float
            (Optional) Target conduction band offset.

        Returns
        -------
        OptmizeResult
            if optimization was succesfull
        float
            -1 is returned if number of layres > 2 (no compensation is performed)
        """
        materials = []
        names = []
        strains = []
        Lmat = {}
        for l in self.layers:
            if l.material.name not in names:
                materials.append(l.material)
                names.append(l.material.name)
                strains.append(l.material.h)
                Lmat[l.material.name] = 0.

            for m in materials:
                if l.material.name == m.name:
                    Lmat[m.name] += l.width

        if len(materials) == 2:

            mat1 = materials[0]
            mat2 = materials[1]

            if x is not None: mat1.updateAlloy(x)
            if y is not None: mat2.updateAlloy(y)

            if Ec is not None:
                DEc = Ec
            else:
                DEc = mat2.params['Ec'] - mat1.params['Ec']
            D1 = Lmat[mat1.name]
            D2 = Lmat[mat2.name]

            minfunc = lambda x : np.sum(np.abs(np.array([DEc,-D2/D1]) - \
                np.array(self.__strainf(x[0],x[1],mat1,mat2))) )

            res = opt.minimize(minfunc, x0=[mat1.x,mat2.x],
                               bounds=[(0,1),(0,1)],
                               method="SLSQP")

            for l in self.layers:
                for im in range(2):
                    if l.material.name == names[im]:
                        l.material.updateAlloy(res.x[im])
            return res

        else:
            print(f"WARNING: Number of materials is {len(materials)}>2!")
            print(f"\tNot implemented: Skipping compensation of strain.")
            return -1

    @staticmethod
    def __strainf(x,y,well,barr):
        """Function for minimizing strain

        Parameters
        ----------
        x : float
            Alloy composition of well material.
        y : float
            Alloy composition of barrier material.
        well : Material
            Well material.
        barr : Material
            Barrier material

        Returns
        -------
        DEc : float
            Conduction band offset for given x and y.
        hwbyhb: float
            h_well/h_barr, where h_i is the substrate lattice mismatch of
            material i.
        """

        x0 = well.x
        y0 = barr.x
        well.updateAlloy(x)
        barr.updateAlloy(y)
        DEc = barr.params['Ec'] - well.params['Ec']
        hwbyhb = well.h/barr.h
        well.updateAlloy(x0)
        barr.updateAlloy(y0)
        return DEc, hwbyhb

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

    def save(self, path = ".", name = None):
        """
        Saves this Structure

        Parameters:
        path : string
            The path to the directory where the material is saved.
        name : string
            The name of this structure to be used when saving

        Returns:
            The file name string if the saved file.
        """
        if name is None:
            name = self.name
            #name = str(self.sid)

        root = etree.Element("Structure")

        author = "[aftershoq] " + str(datetime.datetime.now()) + ""
        content = "This is an autogenerated file."
        header = etree.SubElement(root, "Header")
        etree.SubElement(header, "Author").text = author
        etree.SubElement(header, "Content").text = content
        # Basic information
        etree.SubElement(root, "Name").text = self.name
        etree.SubElement(root, "Nlayers").text = str( len( self.layers ) )
        etree.SubElement(root, "Ndoping").text = str( len( self.dopings ) )
        matlist = []
        layers = etree.SubElement(root, "Layers")
        for l in self.layers:
            # TODO only save if material occures for first time
            filename = l.material.save(path)
            layer = etree.SubElement(layers, "Layer")
            etree.SubElement(layer, "w").text = str(l.width)
            etree.SubElement(layer, "material").text = filename
            etree.SubElement(layer, "lambda").text = str( l.lam )
            etree.SubElement(layer, "eta").text = str( l.eta )
        dopings = etree.SubElement(root, "Doping")
        for d in self.dopings:
            region = etree.SubElement(dopings, "Region")
            etree.SubElement(region, "zi").text = str( d[0] )
            etree.SubElement(region, "zf").text = str( d[1] )
            etree.SubElement(region, "dens").text = str( d[2] )

        filename = name + self.postfix
        file = Path.join(path, filename)
        etree.ElementTree(root).write(file, xml_declaration=True, encoding="utf-8", pretty_print=True)

        return filename


    @staticmethod
    def load(filename):
        """
        Static method for loading a Structure from a file

        Parameters:
        filename: string
            the name of the file to load

        Returns:
            A new Structure
        """

        s = Structure()

        root = etree.parse(filename)
        struct = root.getroot()

        for el in struct:
            if el.tag == "Layers":
                for l in el:
                    for key in l:
                        if key.tag == "w":
                            width = float(key.text)
                        elif key.tag == "material":
                            material = Material.load(key.text)
                        elif key.tag == "eta":
                            eta = float(key.text)
                        elif key.tag == "lambda":
                            lam = float(key.text)
                    s.addLayerIFR(width, material, eta, lam)
            elif el.tag == "Doping":
                for l in el:
                    for key in l:
                        if key.tag == "zi":
                            zi = float(key.text)
                        elif key.tag == "zf":
                            zf = float(key.text)
                        elif key.tag == "dens":
                            dens = float(key.text)
                    s.addDoping(zi, zf, dens)
            elif el.tag == "Name":
                s.name = el.text

        return s

    def get_param(self, param, npoints=1000, nperiods=1):
        """
        Get z and parameter for n periods

        Parameters:
        param : string
            Key to extract parameter
        npoints : int
            Discretization points per period
        nperiods : int
            Number of repetitions of the period

        Returns: Array (z, param)
        """
        zspan = np.linspace(0,self.length,num=npoints,endpoint=False)
        p = []
        for z in zspan:
                ind = self.layerIndex(z)
                p.append(self.layers[ind].material.params[param])

        if nperiods > 1:
            pper = []
            for i in range(nperiods):
                pper = pper + p

        else:
            pper = p

        zspan = np.linspace(0,self.length*nperiods,npoints*nperiods,endpoint=False)

        return np.asarray([zspan,pper])

    def get_conduction_band(self, npoints=100, nperiods=1):
        """
        Get the z and the conduction band for n periods

        Parameters:
        npoints : int
            Discretization points per period
        nperiods : int
            Number of repetitions of the period

        Returns: Array (z, conduction band)
        """
        zspan = np.linspace(0,self.length,num=npoints,endpoint=False)
        cbo = []
        for z in zspan:
                ind = self.layerIndex(z)
                cbo.append(self.layers[ind].material.params["Ec"])

        if nperiods > 1:
            ncbo = []
            for i in range(nperiods):
                ncbo = ncbo + cbo

        else:
            ncbo = cbo

        zspan = np.linspace(0,self.length*nperiods,npoints*nperiods,endpoint=False)

        return np.asarray([zspan,ncbo])



class Material(object):
    '''Defines a material or alloy between two materials.'''

    postfix = ".aqmt"

    # parameter dictionary:
    params_dict = {
        "meff": 0,    # effective mass
        "Ec" : 0,     # Conduction band offset (eV)
        "Eg" : 0,     # Gamma valley band gap (eV)
        "EX" : 0,     # X valley band gap (eV)
        "EL" : 0,     # L valley band gap (eV)
        "Eso": 0,     # Split-off hole band gap (eV)
        "EDel" : 0,   # Delta valley band gap (eV)
        "Ep" : 0,     # Kane energy (eV)
        "del0" : 0,   # Spin-orbit splitting (eV)
        "Valloy" : 0, # Alloy scattering potential (eV)
        "ELO" : 0,    # Longitudinal optical phonon energy (eV)
        "ETO" : 0,    # Transversal optical phonon energy (eV)
        "eps0" : 0,   # Static dielectric constant
        "epsinf" : 0, # High-frequency dielectric function
        "ac" : 0,     # Conduction band deformation potential (eV)
        "acDel" : 0,  # Indirect Delta cond. band def. pot. (eV)
        "acL" : 0,    # Indirect L cond. band def. pot. (eV)
        "av" : 0,   # Valence band deformation potential (eV)
        "c11"  : 0,   # Elastic constant
        "c12"  : 0,   # Elastic constant
        "c44"  : 0,   # Elastic constant
        "vlong" : 0,  # Longitudinal sound velocity (m/s)
        "massdens" : 0,   # Mass density (kg/m^3)
        "molV" : 0,  # Mol volume (nm^3/mol) = alc**3/4 for Zinc Blende
        "lattconst" : 0  # Lattice constant (A)
    }

    def __init__(self,name,params_in=None,mat1=None,mat2=None,C=None,x=None, subs=None):
        '''Constructor. Parameters:
        name:  The name of this material
        params: Material parameters dictionary.
        mat1:   Optional, Material one in new alloy
        mat2:   Optional, Material two in new alloy
        C:      Optional, Bowing parameters for new alloy. All default to 0.
        x:      Optional, relative compositon of mat1 to mat2
        subs:   Optional, substrate material, for strain calculation. If None, then
                substrate is the material itself (no strain)
        '''
        self.name = name
        self.hc = None
        self.substrate = subs
        self.strained = False
        if (params_in is None):
            self.params = Material.params_dict.copy()
        else:
            self.params = params_in.copy()
        self.mat1 = mat1
        self.mat2 = mat2
        if (C is None):
            self.C = Material.params_dict.copy()
        else:
            self.C = C.copy()
        self.x = x
        if (x is not None):
            if x > 1:
                print("ERROR: x > 1 in material creation! Stopping.")
                exit(1)
            self.updateAlloy(x)

    def updateAlloy(self, x, reset_strain = False):
        '''Updates the alloy composition of this alloy material with x
        as the new composition.
        Optional argument reset_strain can be set to True to undo any
        previous strain calculations.
        '''

        if reset_strain:
            self.strained = False

        if(self.x is not None):
            self.hc = None
            self.x = x
            self.params = self.alloy(self.mat1,self.mat2,self.C,self.x)

            if self.strained == True:
                self.strained = False
                self.calcStrain()


    def calcStrain(self):
        '''Calculates the strain effects on band parameters. Corrects
        Eg and Ec accorcing to model-solid theory by [van de Walle PRB 1989]
        and [Gresch Thesis ETH Zuerich 2009].
        '''

        if self.substrate is None:
            # If no substrate is given, then material is assumed relaxed
            return

        # relaxed lattice constant
        ar = self.params["lattconst"]
        # a-parallel. Here, we assume that the material is fully strained
        aII = self.substrate.params["lattconst"]
        # lattice mis-match:
        self.h = (ar-aII)/aII
        # Poisson ratio (assuming substrate is in 001 direction
        D001 = 2*self.params["c12"]/self.params["c11"]
        # epsilon-parallel
        epsII = aII/ar - 1.
        # a-perpendicular
        aL = ar*(1. - D001*epsII)
        # epsilon-perpendicular
        epsL = aL/ar - 1.

        # DeltaOmega/Omega
        DOm = 2*epsII + epsL

        # Valence band shift:
        DEv = self.params["av"]*DOm
        # Conduction band shift (naive):
        if self.strained == False:
            self.params["Ec"] += self.params["ac"]*DOm
        self.strained = True

    def hcrit(self, Nself = 10):
        """
        Returns the critical thickness in Å (calculates it the first time).
        Uses the Matthews and Blakeslee 1974 formula.

        Nself = 10: number of self-consistent iterations. The default is 10, but
        5 also normally gives sufficient precision for evaluation of structures.
        """
        if self.substrate is self or self.substrate is None:
            return 10000000.0

        if self.hc is None:

            b = self.params["lattconst"]/2. # Burgers vector

            hc0 = self.__hrec__(b)
            for i in range(Nself):
                hc0 = self.__hrec__(hc0)

            self.hc = hc0

        return self.hc

    def __hrec__(self, hc):
        """
        Uses Matthews and Blakeslee 1974 for threading dislocation.
        """
        a = self.params["lattconst"]
        asub = self.substrate.params["lattconst"]
        b = a/np.sqrt(2.) # Burgers vector
        f = np.abs(a-asub)/asub # Lattice mismatch
        # Poisson ratio:
        nu = self.params["c12"]/(self.params["c11"]+self.params["c12"])

        return b*(1-nu/4)/(4*np.pi*f*(1+nu))*(np.log(hc/b)+1)


    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)

    def save(self, path = ".", name = None):
        """
        Saves this Material

        Parameters:
        path : string
            The path to the directory where the material is saved.

        Returns:
            The file name string if the saved file.
        """

        if name is None:
            name = self.name

        root = etree.Element("Material")

        # Header:
        author = "[aftershoq] " + str(datetime.datetime.now()) + ""
        content = "This is an autogenerated file."
        header = etree.SubElement(root, "Header")
        etree.SubElement(header, "Author").text = author
        etree.SubElement(header, "Content").text = content
        # Basic information
        etree.SubElement(root, "Name").text = self.name
        etree.SubElement(root, "Strained").text = str(self.strained)
        etree.SubElement(root, "Alloy").text = str(self.mat1 is not None)
        fileSubst = "None"
        if self.substrate is not None:
            fileSubst = self.substrate.save(path)
        etree.SubElement(root, "Substrate").text = fileSubst

        parameters = etree.SubElement(root, "Params")
        for key in self.params:
            etree.SubElement(parameters, key).text = str(self.params[key])

        bowings = etree.SubElement(root, "Bowing")
        for key in self.C:
            etree.SubElement(bowings, key).text = str(self.C[key])

        if self.mat1 is not None:
            filemat1 = self.mat1.save(path)
            etree.SubElement(root, "Mat1").text = filemat1

        if self.mat2 is not None:
            filemat2 = self.mat2.save(path)
            etree.SubElement(root, "Mat2").text = filemat2

        if self.x is not None:
            etree.SubElement(root, "x").text = str(self.x)


        # Print to file
        filename = name + self.postfix
        file = Path.join(path, filename)
        etree.ElementTree(root).write(file, xml_declaration=True, encoding="utf-8", pretty_print=True)

        return filename

    @staticmethod
    def load(filename):
        """
        Static method for loading a Material from a file

        Parameters:
        filename: string
            the name of the file to load

        Returns:
            A new Material instance
        """
        root = etree.parse(filename)
        mat = root.getroot()
        params = {}
        C = {}
        mat1, mat2 = None, None
        x = None
        for el in mat:
            if el.tag == "Name":
                name = el.text
            elif el.tag == "Substrate":
                if el.text != "None":
                    substrate = Material.load(el.text)
                else:
                    substrate = None
            if el.tag == "Params":
                parEl = el
            elif el.tag == "Bowing":
                parBow = el
            elif el.tag == "Alloy":
                alloy = el.text
            elif el.tag == "Mat1":
                if el.text != "None":
                    mat1 = Material.load(el.text)
            elif el.tag == "Mat2":
                if el.text != "None":
                    mat2 = Material.load(el.text)
            elif el.tag == "x":
                x = float(el.text)

        for el in parEl:
            try:
                val = float(el.text)
            except ValueError:
                pass
            params.update({el.tag : val})

        for el in parBow:
            try:
                val = float(el.text)
            except ValueError:
                pass
            C.update({el.tag : val})

        return Material(name,params_in=params,mat1=mat1,mat2=mat2,C=C,x=x, subs=substrate)

    @staticmethod
    def alloy(mat1,mat2,C,x):
        '''Creates and returns a new set of material parameters, which are
        defined by the parameters of mat1, mat2, the bowing parameters
        in C, and the composition x. C is a list of bowing parameters matching
        the parameters in Material.params and MaterialPar:
        newparams = x*param1 + (1-x)*param2 - C*x*(1-x)
        x is the alloy fraction of mat1 to mat2.
        '''

        params3=Material.params_dict.copy()
        for key in mat1.params:
            params3[key]=(mat1.params[key]*x + mat2.params[key]*(1-x) -
                          C[key]*x*(1-x))
        # Alloy scattering:
        params3["Valloy"] = x*(1-x)*(mat1.params["Ec"]-mat2.params["Ec"])
        return params3
