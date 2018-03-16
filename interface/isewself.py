'''
Created on 12 Mar 2018

@author: martin
'''

from interface import Interface
from classes import MaterialPar as mp
import const

class Isewself(Interface):
    '''
    Interface for sewself.
    '''
    # dictionary with name and default value
    sewselfpar = {
        "emin":-0.0001,
        "emax":0.14,
        "n":1200,
        "stepwf":3,
        "iniz0":30,
        "fin0":25,
        "writecoef":1,
        "gnuplotgenere":1,
        "eerr0":1e-8,
        "very_small":1e-21,
        "incmin":1e-20,
        "clockmax":100,
        "clockmin":7,
        "Eabove":0.14,
        "rightmost":-100,
        "leftmost":-400,
        "mumin":-0.17,
        "mumax":0.17,
        "stepmu":0.002,
        "imax":12,
        "rho":0.3,
        "itarget":1,
        "idepth":3,
        "convergence":1e-7,
        "Espan":0.3,
        "maxiter":50,
        "HWHM":0.0001,
        "nref":3.6,
        "lambda":0.080
        }
    
    matpar = {
        "efield" : 0.0,
        "ldiff"  : 0.0,
        "nbswellinj":4,
        "nbswellact":4,
        "coeffls":1,
        "bool_inj":0,
        "bool_act":0,
        "bool_one":1,
        "matchoice":1,
        "inc0":5,
        "incw":2
        }

    def __init__(self, binpath, pltfm, numpar, material_list):
        '''
        Constructor
        '''
        super(Isewself,self).__init__(binpath, pltfm, numpar)
        self.material_list = material_list
        
    def runStructures(self,structures,numpar,path):
        for ss in structures:
            self.initdir(ss)
            
    def initdir(self,structure,path):
        spath = path + str(structure.sid)
    
    def gatherResults(self,structures,path):
        pass
    
    def writeStructFile(self,structure,path):
        filepath = path + "/structure.txt"
        with open(filepath,'w') as f:
            f.write(str(self.matpar["efield"]) + " # efield\n")
            f.write(str(self.matpar["ldiff"]) + " # diffusion length\n")
            f.write(str(len(structure.layers)) + " # number of layers\n")
            f.write(str(self.matpar["nbswellinj"]) + " # nbswellinjector\n")
            f.write(str(self.matpar["nbswellact"]) + " # nbswellactive\n")
            f.write(str(self.matpar["coeffls"]) + " # coeff (longer/shorter structure)\n")
            f.write("switches for structure length\n")
            f.write(str(self.matpar["bool_inj"]) + " # injector only\n")
            f.write(str(self.matpar["bool_act"]) + " # active only\n")
            f.write(str(self.matpar["bool_one"]) + " # one period\n")
            f.write("material parameters\n")
            f.write(str(self.matpar["matchoice"]) + " matchoice\n")
            for mat in self.material_list:
                f.write(str(mat) + " ")
                f.write(str(mat.params[mp.Ec]) + " ")
                if mat.x is None:
                    f.write(str(1))
                else:
                    f.write(str(mat.x))
                f.write("\n")
            f.write( "end discretisation of the potential\n" )
            f.write( str(self.matpar["inc0"]) + " inc0\n" )
            f.write( str(self.matpar["incw"]) + " incw\n" )
            f.write( "Structure (mat, thick (A), doping x 1e-18) \n" )
            i = 0
            for layer in structure.layers:
                f.write(str(layer.material) + " " + str(layer.width) + " " )
                doping = 0
                for dop in structure.dopings:
                    if dop[0]>=structure.layerPos(i):
                        doping += dop[2]
                f.write( str( doping*1e-18 ) + "\n")


    
    def writeParameterFile(self,path):
        filepath = path + "/mat.par"
        with open(filepath,'w') as f:
            f.write("mwell, mbarrier, Gamma, hlo, kp0\n")
            wellind = 0
            for i in range(1,len(self.material_list)):
                if self.material_list[i].params[mp.Ec] < self.material_list[wellind].params[mp.Ec]:
                    wellind = i
            if wellind < len(self.material_list)-1:
                barrind = wellind + 1
            else:
                barrind = wellind -1
            well = self.material_list[wellind]
            barr = self.material_list[barrind]
            f.write(str(well.params[mp.meff]) + " ")
            f.write(str(barr.params[mp.meff]) + " ")
            gamma = const.hbar_eV**2*const.qe/(2*const.me)
            gamma = gamma / well.params[mp.Eg] / well.params[mp.meff]
            f.write(str(gamma) + " ")
            f.write(str(well.params[mp.ELO]) + " ")
            epsp = 1./well.params[mp.epsinf] - 1./well.params[mp.eps0]
            f.write(str(1./epsp) + " ")
            f.write(str(well) + "/" + str(barr) + "\n")
            
    def writeSewselfPar(self,path):
        filename = path + "/sewself.par"
        with open(filename, 'w') as f:
            d = self.sewselfpar
            order = [
                "emin",
                "emax",
                "n",
                "stepwf",
                "iniz0",
                "fin0",
                "writecoef",
                "gnuplotgenere",
                "eerr0",
                "very_small",
                "incmin",
                "clockmax",
                "clockmin",
                "Eabove",
                "rightmost",
                "leftmost",
                "mumin",
                "mumax",
                "stepmu",
                "imax",
                "rho",
                "itarget",
                "idepth",
                "convergence",
                "Espan",
                "maxiter",
                "HWHM",
                "nref",
                "lambda"
                ]
            for key in order:
                value = d.get(key)
                f.write(str(value) + " " + str(key) + "\n")
                        