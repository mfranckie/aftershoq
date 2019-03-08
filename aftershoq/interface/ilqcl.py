'''
Created on 16 Mar 2018

@author: Martin Franckie

'''

from aftershoq.interface import Interface
from aftershoq.structure import Structure
import aftershoq.utils.systemutil as su
from aftershoq.numerics.runplatf import Local
import time
import aftershoq.utils.debug as dbg
import subprocess
from matplotlib import pyplot as pl
import numpy as np
from aftershoq.utils import const
from aftershoq.materials import GaAs
from path import Path
import json

class Ilqcl(Interface):
    '''
    Inteface for the LQCL program package, documented in:
    [KirsanskasPRB2018]
    '''

    proglqcl = "run_lqcl.py"


    def __init__(self,binpath="",pltfm=Local(),wellmaterial=GaAs()):
        '''
        Constructor. Subclass specific parameters:
        binpath: Path to the calcWS program executable
        wellmaterial: Material of the well of the structure (used for
        dielectric and phonon properties).
        '''

        # Input to calcWS program
        negf_numpar = {
            "Npint" : 5,
            "cshift" : 0,
            "Lbound" : 1.E-3,
            "Ubound" : 1.E-4,
            "Blochtol" :1.E-6,
            "Nh" : 0,
            "Igauge": 1,
            "Eminsub" : 0,
            "Emaxadd": 0,
            "gen" : 1.E-4,
            "Iconv": 2,
            "Bei":0.8,
            "Nhist": 40,
            "boolPrinc" : False,
            "boolEins" : False
            }

        # Parameters for lqcl:
        lqcl_numpar = {
            "neighbours_imp" : 1,
            "neighbours_ifr" : 2,
            "neighbours_imp" : 1,
            "neighbours_z" : 1,
            "neighbours_LO" : 1,
            "use_drive_field" : False,
            "use_wannier_dir" : False,
            "keep_wsdir" : True,
            "qmin" : 1e-10,
            "qmax" : 3.,
            "qzmin" : -3.,
            "qzmax" : 3.,
            "prntq" : False,
            "zshift" : 0.
        }

        super(Ilqcl,self).__init__(binpath, pltfm)
        self.numpar.update(negf_numpar)
        self.numpar.update(lqcl_numpar)
        self.progWS = Path.joinpath(binpath,"calcWS.out")


        self.wellmat = wellmaterial
        self.processes= []
        self.datpath = "data"
        self.wannierpath = "wannier"
        self.numparfile = "numpar_lqcl.json"

    def __str__(self):
        return "Ilqcl"

    def initdir(self,ss,path):
        '''Initialize the dierctory for Structure s, with base path "path".'''
        pathdata=Path.joinpath(path,self.datpath)
        pathwannier=Path.joinpath(path,self.wannierpath)
        su.mkdir(pathdata)
        su.mkdir(pathwannier)
        self.writeWannier(ss,pathwannier)
        self.writeMaterial(self.wellmat, "# "+str(self.wellmat.name),pathwannier)
        self.writeNumpar(path)
        self.wellmat.save(path, "wellmaterial")
        ss.save(path)

    def writeNumpar(self, path = "."):
        """
        Writes self.numpar into file named self.numparfile

        Attributes:
        path : string
            directory in which to write.
        """

        f = open(Path.joinpath( path, self.numparfile ) , "w" )
        # for key in self.numpar:
        #     f.write(key + " " + str(self.numpar[key]) + "\n")
        json.dump(self.numpar, f, indent = 4)
        f.close()




    def runStructures(self,structures,path, calcStates = True):
        '''Run simulations for all structures in the given structure list with
        the base path "path". This method dispatches all processes and returns
        the user has to wait for processes to finish before accessing results.

        Stores started processes in self.processes
        '''

        local = Local()
        for ss in structures:
            spath = Path.joinpath( path, str(ss.dirname) )
            su.mkdir(spath)
            if calcStates:
                self.initdir(ss, spath)
                #proc = su.dispatch(self.progWS, "",spath)
                proc = self.pltfm.submitjob(self.progWS,[],spath,1,"00:10")
                self.processes.append(proc)
        if calcStates:
            dbg.debug("Starting calcWS program.....\n",dbg.verb_modes["verbose"],self)
            dbg.flush()
        # TODO do not wait for all processes-start each strucrure when ready!
        self.waitforproc(1)
        dbg.debug("Starting lqcl.....\n",dbg.verb_modes["verbose"],self)
        dbg.flush()
        #del processes
        #processes = []
        for ss in structures:
            spath = Path.joinpath( path , str(ss.dirname) )

            proc = self.pltfm.submitjob(self.proglqcl,[],Path.joinpath(spath,self.datpath) )

            self.processes.append(proc)
        return self.processes



    def checkactive(self):
        pactive = False
        for p in self.processes:
                if self.pltfm.jobstatus(p):
                    pactive=True
        return pactive


    def gatherResults(self, structures, pathwd, pathresults = None, runprog = True):
        '''Write results to pathresults/results.log and run hdiag and bandplot
        in pathwd/s.dirname/self.datpath/eins/x/ for each i and x. Stores WS
        resutls as a new attribute levels[directory][WS level][data field] in
        each Structure object in the list structures.
        '''

        if(pathresults is None):
            pathresults = pathwd

        with open(pathresults+'/results.log','a') as f:
            f.write('# Results for structures:\nID | merit | \
            N times layer width | Ndop times (zi, zf, nvol) | N times x \n')
            for ss in structures:
                spath = Path.joinpath( pathwd , str(ss.dirname) )

                ss.merit = self.getMerit(ss, pathwd)

                f.write(str(ss.dirname) + " ")
                f.write(str(ss.merit) + " " )
                for layer in ss.layers:
                    f.write(str(layer.width)+" ")
                for doping in ss.dopings:
                    for val in doping:
                        f.write( str(val) +" ")
                for layer in ss.layers:
                    f.write(str(layer.material.x)+" ")

                f.write("\n")
                if self.merit == self.merits["Chi2"]:
                    su.mkdir(pathresults + ss.dirname)
                    with open(pathresults + ss.dirname+"/chi2.log", 'w') as chif:
                            chif.write("Not implemented yet in Ilqcl!")

    def getMerit(self,structure,path):
        '''Returns the merit function evaluated for the Structure structure,
        with base path "path".
        '''

        # First try to see if results were already evaluated:
        # assumes structure for results.log:
        # sid | merit | ...

        if 'results.log' in su.ls(path):
            with open(path + '/results.log') as f:
                found = False
                for line in f:
                    line = line.split()
                    if line[0] == str(structure.sid):
                        return line[1]

        path = path + "/" + structure.dirname + self.datpath

        results = []
        out = 0.

        if self.merit==Interface.merits.get("max gain") :
            pass
        elif self.merit == Interface.merits.get("(max gain)/(current density)"):
            pass
        elif self.merit == self.merits.get("estimated gain") :
            pass
        elif self.merit == self.merits["Chi2"]:

            if len( self.target ) > 0:
                E1 = self.target[0]
            else:
                E1  = None
            if len(self.target ) > 1:
                E2 = self.target[1]
            else:
                E2 = None

            if len( self.target ) > 2:
                gamma = self.target[2]
            else:
                gamma = None

            chi2 = self.calcChi2(structure, E1, E2, gamma)
            if chi2 < 0:
                return "ERROR"
            else:
                return chi2
        elif self.merit == self.merits["photocurrent"]:
            omegat = self.target[0]
            if len(self.target) > 1:
                a = self.target[1]
            else:
                a = 0.010 # eV
        else:
            print("No such merit function!")

        return out

    def writeMaterial(self,material,nametag,dirpath = None):
        '''Writes the material.inp input file.'''

        if dirpath is None:
            path="material.inp"
        else:
            path = dirpath + "/material.inp"
        with open(path, 'w') as f:
            f.write(nametag + "\n")
            f.write(str(material.params["meff"])+" # meff\n")
            f.write(str(material.params["eps0"])+" # eps0\n")
            f.write(str(material.params["epsinf"])+" # epsinf\n")
            f.write(str(material.params["ELO"])+" # ELO\n")
            f.write(str(material.params["ac"])+" # deformation potential\n")
            f.write(str(material.params["vlong"])+" # vlong\n")
            f.write(str(material.params["massdens"])+" # mass density\n")
            f.write(str(material.params["molV"])+ " # mol volume")
        f.closed

    def writeWannier(self,struct,dirpath=None):
        '''Writes the wannier8.inp input file.'''

        if dirpath is None:
            path="wannier8.inp"
        else:
            path = dirpath + "/wannier8.inp"

        # re-scaling of the CBO to lowest energy:
        Ec = []
        [Ec.append(l.material.params["Ec"]) for l in struct.layers]
        Ecmin = min(Ec)


        try:
            with open(path,'w') as f:
                f.write(str(struct.Nl))
                f.write("\n# NS lines follow below containing 4 parameters for each layer")
                f.write("\n# DS: width of layer (nm)\n# EcS: Conduction band at Gamma point (eV),")
                f.write("\n# MS: effective mass at Gamma point in units of me")
                f.write("\n# AlloyS: =x*(1-x)*[Ec(A)-Ec(B)] (eV) for ternary lattices A_xB_(1-x)C\n")
                il = 1
                for l in struct.layers:
                    f.write(str(l.width)+" ")
                    f.write(str(l.material.params["Ec"]-Ecmin)+ " ")
                    f.write(str(l.material.params["meff"])+ " ")
                    f.write(str(l.material.params["Valloy"])+ " ")
                    f.write(" # Layer " + str(il))
                    il=il+1
                    f.write("\n")
                f.write(str(self.wellmat.params["Ep"])+"\n")
                f.write(str(len(struct.dopings))+"\n")
                f.write("# Left and right end point of doped region (nm)\n"\
                        "# Doping density per volume (1/cm^3)\n"\
                        "# Number of scattering layers for each region\n")
                il = 1
                for d in struct.dopings:
                    f.write(str(d[0])+" "+str(d[1])+" "+"{0:E}".format(d[2])+" "+ str(1)+" # region " + str(il) +"\n")
                    il = il+1
                f.write(str(struct.layers[0].eta) + " " + str(struct.layers[0].lam) + " # eta, lambda\n")
                f.write(str(self.numpar["Nstates"]) + "\t\t# Number of Wannier states per period to be calculated\n")
                f.write(str(self.numpar["Nz"]) + "\t\t# Number of z points per period for wave functions\n")
                f.write(str(self.numpar["Nq"]) + "\t\t# Number of q points\n")
                f.write(str(self.numpar["Nper"])+ "\t\t# Number of periods for wave functions in output\n")
                f.write("# The following parameters may treat numerical problems - usually unchanged\n")
                f.write(str(self.numpar["Npint"])+"\t\t# increases the internally used number of periods\n")
                f.write(str(self.numpar["cshift"])+"\t\t# shift between symmetry point and center of period\n")
                f.write(str(self.numpar["Lbound"]) +" " + str(self.numpar["Ubound"])+ "\t\t# lower bound for energies/gaps of minibands\n")
                f.write(str(self.numpar["Blochtol"])+ "\t\t# Tolerance for accuracy of solution for Blochfunction\n")
            f.closed
        except IOError:
            print("WARNING: Directory "+dirpath+" not found!")

    def loadStructures(self, origs, path):
        '''
        Load structures from directory tree. Assumes the following structure:
        path/sid/wannier8.inp/, where sid is the structure id (integer).

        Parameters

        origs: model structure, where layer widths will be replaced by read
        widths, for each structure
        path: path to root directory of tree.
        '''

        structures = []
        dirs = su.listdirs(path)
        for folder in dirs:
            s = Structure(origs)
            with open(path+"/"+folder+"/wannier8.inp", 'r') as wannier:
                s.sid = int( folder )
                s.dirname = str( s.sid )
                Nl = int( wannier.readline() )
                Ncomm = 5
                for _ in range(Ncomm):
                    wannier.readline()
                for i in range(Nl):
                    lw = float( wannier.readline().split()[0] )
                    s.layers[i].width = lw
            structures.append(s)
        return structures



    def plotbands(self, structure, path, einspath = None):
        '''
        Plot the band structure and Wannier-Stark states for a structure.
        Must first have run gatherResults() for structure.

        Parameters

        structure: Structure for which to plot the bands.
        path: The path to the base of the tree
        einspath (optional): specifies which parameters to plot for
        '''

        pass

    def calcChi2(self, structure, E1, E2, gamma = None):
        """
        Calcualtes the second order susceptibility Chi(2) for the three-level
        system |0>, |1>, |2>, where E_2 > E_1 > E_0 according to
        Dupont et al, IEEE J. Quant. Electron. 42, pp. 1157-1174 (2006)

        Parameters:
        structure: The structure to calculate Chi(2) for
        E1: The pump photon energy hbar*omega (eV). If E1 is None, then
        this is taken to be E_2-E_0
        E2: The idler photon energy hbar*omega (eV). If E2 is None, then
        this is taken to be E_1-E_0
        gamma: The FWHM broadening to be used in the calculation for Chi(2).
        If gamma is None, then this is calculated from the level broadenings.

        """
        pass
