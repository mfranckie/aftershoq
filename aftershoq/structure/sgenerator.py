'''Created on 29 Mar 2018

@author: Martin Franckie
'''


from aftershoq.structure import Structure

import random
from hilbert_curve.hilbert import HilbertCurve
from aftershoq.utils import HilbertUtil
import numpy as np

doping_modes = {"FOLLOW LAYER SHEET": 0,
                "FOLLOW LAYER VOLUME": 1,
                "FIXED POSITION SHEET": 2,
                "FIXED POSITION VOLUME": 3}

comp_modes = {"UNLINKED":0,
              "LINKED LAYERS": 1,
              "LINKED ALL": 2}

width_modes = {"CONTINUOUS": 0,
               "MONOLAYER" : 1,
               "DISCRETE"  : 2}

strain_modes = {"NO STRAIN" : 0,
                "COMPENSATE X" : 1}

class Sgenerator():

    def __init__(self,origstruct,dw=None,dx=None,ddop=None,
                 doping_mode = doping_modes.get("FOLLOW LAYER SHEET"),
                 comp_mode = comp_modes.get("UNLINKED"),
                 width_mode = width_modes.get("CONTINUOUS"),
                 strain_mode = strain_modes.get("NO STRAIN"),
                 p = 3 ):
        '''
        Constructor which takes:
        origstruct: original structure to modify, with N layers
        dw: list of length N with one range of widths for each layer
            The widths will vary in the interval w0-dw to w0+dw
        dx: list of length N with one range of composition for each layer
            The comp. will vary in the interval x0-dx to x0+dx
            If comp_mode is specified to "LINKED LAYERS", in the end of the
            list must follow M (number of dx[i] > 0) integers with tags
            specifying which layers are to be linked together
            (their x will always change together).
        ddop: (optional) a list of the structure [dpos, dw, ddens]
            where dops gives the range of positions, dw the range of
            doping layer widths, and ddens the range of doping densities.
        doping_mode: (optional) Tells generator how to treat doping layers.
            Can be the following (doping_modes):
            FOLLOW LAYER SHEET*: Doping density follows original sheet dop.
            FOLLOW LAYER VOLUME: Doping density follows original volume dop.
            FIXED POSITION SHEET: Dop. remains at a fixed postion and sheet dop.
            FIXED POSITION SHEET: Dop. remains at a fixed postion and vol. dop.
        comp_mode: (optional) Tells generator how to treat composition changes.
            Can be the following (comp_modes):
            UNLINKED*: Each specified layer's x changes independently.
            LINKED LAYERS: The specified layers' x changes synchronously.
            LINKED ALL: All specified layers' x changes synchronously.
            In the laset two cases, only the ranges of the first element of dx
            is used.
        width_mode: (optional) Tells generator to use continuous or discretized
            widths when generating structures.
            Can be the following (width_modes):
            CONTINUOUS*: Width can vary continuous
            MONOLAYER: Width can vary in monolayers only.
        strain_mode: (optional) Tells generator which strain mode to use:
            NO STRAIN*: No strain compensation is done
            COMPENSATE X: The strain is compensated by changing the alloy
            composition of the layers. Attempts to keep band offset the same.
        p: (optional) Parameter for Hilbert Curve density. Needed if
            HilbertCurve will be used to generate structures.
        '''

        self.doping_mode = doping_mode
        self.comp_mode = comp_mode
        self.width_mode = width_mode
        self.strain_mode = strain_mode

        if dx is None:
            dx = []
            [dx.append(0.) for _ in range(origstruct.Nl)]

        if dw is None:
            dw = []
            [dw.append(0.) for _ in range(origstruct.Nl)]

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


        # Initialize the varying fields:

        # Number of dimensions in Euclidian space
        ND = 0

        # Define list of changeable parameters
        params = []
        dparams = []
        windex = []
        xindex = []
        dopindex = []

        # 1) Layer widths:
        for i in range(0,len(self.orig.layers)):
            if self.dw[i]>0:
                params.append(self.orig.layers[i].width)
                dparams.append(self.dw[i])
                windex.append(i)
                ND += 1

        # 2) Doping layers:
        for d in self.orig.dopings:
            for i in range(0,len(d)):
                if self.ddop[i] > 0:
                    params.append(d[i])
                    dparams.append(self.ddop[i])
                    dopindex.append(i)
                    ND +=1

        # 3) Alloy composition x:
        if self.comp_mode == comp_modes.get("LINKED LAYERS"):
            tags = []
            self.tagindex = []
            itag = 0
            for i in range(0, len(self.orig.layers)):
                if self.dx[i] > 0:
                    tag = self.dx[len(self.orig.layers)+itag]
                    itag += 1
                    if tag not in tags:
                        tags.append(tag)
                        params.append(self.orig.layers[i].material.x)
                        dparams.append(self.dx[i])
                        ND += 1
                    xindex.append(i)
                    self.tagindex.append(tag)
        else:

            for i in range(0,len(self.dx)):
                if self.dx[i] > 0:
                    params.append(self.orig.layers[i].material.x)
                    dparams.append(self.dx[i])
                    xindex.append(i)
                    ND += 1
                    if self.comp_mode == comp_modes.get("LINKED ALL"):
                        # only include one dx in parameter space
                        for j in range(i,len(self.dx)):
                            if self.dx[j] > 0:
                                xindex.append(j)
                        break


        # 4) Add further parameters at this point...

        # define the Hilbert curve:
        hilbert_curve = HilbertCurve(p, ND)

        # imax = length of hilbert curve (# nodes)
        imax = 2**(ND*p)-1

        # pmax = side of hyper cube
        pmax = 2**p-1

        print('Dim = {0}, p = {1}, imax = {2}, pmax = {3}'.format(ND,p,imax,pmax))

        self.hutil = HilbertUtil(hilbert_curve)
        self.params = params
        self.dparams = dparams
        self.windex = windex
        self.dopindex = dopindex
        self.xindex = xindex

    def genRanStructs(self,N):
        """
        Generate N structures with random variations in layer widths.
        Doping modes "FOLLOW LAYER SHEET/VOLUME" assume that doping
        covers entire layer width.
        Doping modes "FIXED WIDTH SHEET/VOLUME" changes doping pos.
        independently of layer widths.
        TODO: Change compositions!
        """
        for _ in range(0,N):
            news = Structure(self.orig)

            for il in range(len(self.dw)):
                w0 = self.orig.layers[il].width
                width = w0 + 2*(random.random()-0.5)*self.dw[il]
                news.layers[il].width = width

            news.dopings = []
            for dl in self.orig.dopings:

                zc = (dl[0]+dl[1])/2
                layerindex = self.orig.layerIndex(zc)
                sheet = self.orig.layers[layerindex].width*dl[2]


                if (self.doping_mode == doping_modes.get("FOLLOW LAYER SHEET") or
                    self.doping_mode == doping_modes.get("FOLLOW LAYER VOLUME")):

                    zi = 0
                    zf = news.layers[layerindex].width

                    if self.doping_mode == doping_modes.get("FOLLOW LAYER SHEET"):
                        dopdens = sheet/news.layers[layerindex].width
                    else:
                        dopdens = dl[2]

                    dopdens = dopdens + 2*(random.random()-0.5)*self.ddop[2]

                else:

                    dzc = zc - self.orig.layerPos(layerindex)

                    # define new offset from layer start
                    dzc = dzc + 2*(random.random()-0.5)*self.ddop[0]
                    dopw = (dl[1]-dl[0]) + 2*(random.random()-0.5)*self.ddop[1]

                    # add new doping layer
                    zi = dzc - dopw/2
                    zf = dzc + dopw/2

                    if self.doping_mode == doping_modes.get("FIXED POSITION SHEET"):
                        dopdens = sheet/(zf-zi)
                    else:
                        dopdens = dl[2]

                    dopdens = dopdens + 2*(random.random()-0.5)*self.ddop[2]

                news.addDoping(zi, zf, dopdens, layerindex)

            if self.strain_mode == strain_modes.get("COMPENSATE X"):
                news.compensate_strain()

            self.structures.append(news)
            del news

    def genRanHilbertStructs(self,N, p=None):
        """
        Generate N structures on the Hilbert curve defined by p (set in
        constructor).
        Number of dimensions is determined by # parameters to
        change. Includes the end-points of the Hilbert curve.
        Returns the parameters scaled to the hyper cube space.
        """

        if p is not None:
            # override p from conscrtuctor
            self.p = p
            ND = self.hutil.hilbert_curve.n
            self.hutil = HilbertUtil(HilbertCurve(p,ND))
            imax = self.hutil.imax
            pmax = self.hutil.pmax
            print('Dim = {0}, p = {1}, imax = {2}, pmax = {3}'.format(ND,p,imax,pmax))

        d = []

        for i in range(0,N):

            # pick a random point along the hilbert curve, and include
            # end points:
            if i == 0:
                d.append(i)
            elif i == N-1:
                d.append(self.hutil.imax)
            else:
                d.append(random.random()*self.hutil.imax)

        coordinates = self.gen_struct_from_hilbert_curve(d)

        return coordinates

    def gen_struct_from_hilbert_curve(self,newd):
        '''
        Generates new structures at the points along the
        Hilbert curve given in newd, and returns their
        scaled coordinates.
        '''

        dpar = np.array(self.dparams)
        parmin = np.array(self.params) - dpar
        ranges = []
        [ranges.append([self.params[i]-dpar[i],self.params[i]+dpar[i]]) for i in range(len(self.params)) ]
        coordinates = []
        for i in range(0,len(newd)):

            # interpolate to get the coordingates:
            p_ipl = self.hutil.interp_coords_from_dist(newd[i])

            # un-scale it and create structure
            parnew = np.array(p_ipl)

            #par_unscaled = parnew*2*dpar/self.hutil.pmax + parmin
            par_unscaled = np.squeeze(
                self.hutil.scale_coordinates(parnew,dpar,parmin) )

            # generate strucure:
            news = Structure(self.orig)
            pindex = 0

            # Layer widths:
            for i in range(0,len(self.windex)):
                news.layers[self.windex[i]].width = par_unscaled[pindex]
                pindex+=1

            news.dopings = []

            dop_start = 0
            dop_step = int(len(self.dopindex)/len(self.orig.dopings))
            for d in self.orig.dopings:

                # Doping layers:
                z0 = d[0]
                zend = d[1]
                doping = d[2]
                center = (z0 + zend)/2
                width = (zend - z0)
                sheet = width*doping
                zizfdens = [center, width, doping]

                dop_end = dop_start + dop_step

                for i in range(dop_start,dop_end):

                    zizfdens[self.dopindex[i]] = par_unscaled[pindex]

                    pindex +=1

                dop_start = dop_end

                if (self.doping_mode == doping_modes.get("FOLLOW LAYER SHEET") or
                    self.doping_mode == doping_modes.get("FOLLOW LAYER VOLUME")):

                    layerindex = self.orig.layerIndex(center)
                    z0 = news.layerPos(layerindex)
                    zend = z0 + news.layers[layerindex].width

                    if (self.doping_mode == doping_modes.get("FOLLOW LAYER SHEET") ):
                        doping = zizfdens[2]*width/(zend-z0)
                    else:
                        doping = zizfdens[2]

                else:
                    z0 = zizfdens[0]
                    zend = zizfdens[1]

                    if (self.doping_mode == doping_modes.get("FIXED POSITION SHEET") ):
                        doping = zizfdens[2]/width*(zend-z0)
                    else:
                        doping = zizfdens[2]

                    doping = zizfdens[2]

                news.addDoping(z0,zend,doping)

            # Alloy composition:
            tags = []
            ptag = []
            for i in range(0,len(self.xindex)):


                if self.comp_mode == comp_modes.get("UNLINKED"):
                    news.layers[self.xindex[i]].material.updateAlloy(par_unscaled[pindex])
                    pindex +=1
                elif self.comp_mode == comp_modes.get("LINKED ALL"):
                    news.layers[self.xindex[i]].material.updateAlloy(par_unscaled[pindex])
                    pass
                elif self.comp_mode == comp_modes.get("LINKED LAYERS"):
                    print("Layers are linked:")
                    tag = self.tagindex[i]
                    print("tag = ", tag)
                    if tag not in tags:
                        tags.append(tag)
                        ptag.append(pindex)
                    else:
                        pindex = ptag[tags.index(tag)]

                    print("pindex = ", pindex)
                    news.layers[self.xindex[i]].material.updateAlloy(par_unscaled[pindex])
                    news.layers[self.xindex[i]].material.name += str(tag)


                    if len(ptag) > 0:
                        pindex = np.max(ptag) + 1
                    else:
                        pindex += 1

            if self.width_mode == width_modes.get("MONOLAYER"):
                news.convert_to_ML()
            if self.strain_mode == strain_modes.get(["COMPENSATE X"]):
                news.compensate_strain()

            self.structures.append(news)

            coordinates.append(parnew)

        return coordinates


    def load_structs_on_hilbert_curve(self, structures):

        self.structures = structures

        return self.coords_from_struct(structures)


    def coords_from_struct(self, struct):

        if struct.__class__ is not list:
            struct = [struct]


        coords = []

        for s in struct:

            params = []

            # 1) Layer widths:
            for i in self.windex:
                params.append(s.layers[i].width)

            # 2) Doping layers:
            for d in s.dopings:
                for i in self.dopindex:
                    params.append(d[i])

            # 3) Alloy composition x:
            if self.comp_mode == comp_modes.get("LINKED LAYERS"):
                for i in self.xindex:
                    tag = self.dx[len(s.layers)+itag]
                    itag += 1
                    if tag not in tags:
                        tags.append(tag)
                        params.append(s.layers[i].material.x)

            else:
                for i in self.xindex:
                    params.append(s.layers[i].material.x)


            params = np.array(params)
            dpar = np.array(self.dparams)
            parmin = np.array(self.params) - dpar

            coords.append(np.squeeze(self.hutil.unscale_coordinates(params,dpar,parmin)))

        return coords
