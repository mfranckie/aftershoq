'''
Created on 12 Mar 2018

@author: Martin Franckie
'''

from .interface import Interface
import structure.matpar as mp
from utils import const
import utils.systemutil as su
import time
import utils.debug as dbg
import subprocess
import numpy as np

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
    
    numpar = {
        "efield" : 0.0,
        "lattice_temp" : 77,
        "el_temp" : 120,
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

    def __init__(self, binpath, pltfm, material_list, commands = None):
        '''
        Constructor. Baseclass-specific parameters:
        material_list: A list of the Materials that will be used.
        '''
        
        super(Isewself,self).__init__(binpath, pltfm)
        self.material_list = material_list
        self.processes = []
        self.buildmatself = binpath + "buildmatself_MAC"
        self.sewself = binpath + "sewself.2.32_MAC"
        self.structfilename = "structure.txt"
        self.merits.update({'DeltaE_12' : 8, 'Elase' : 9, 'QWIP' : 10})
        self.merits.update({'absorption': 11})
        
        if commands == None:
            commands = ['e', 'd', 'c', 'a210']
#             commands = ""
#             commands  += "e\n"   # compute self-consitent potential
#             commands += str(self.numpar.get('lattice_temp')) + "\n"
#             commands += str(self.numpar.get('el_temp')) + "\n"
#             commands += "d\n"   # compute abosrption from ground state
#             commands += "c\n"   # compute qwip absorption from upper states
#             commands += "a\n2\n1\n0\n" # compute dipole and lifetime between two states
#             commands +=  "q\n"
        
        self.commands = commands
            
        
    def __str__(self):
        return "Sewself"
        
    def runStructures(self, structures, path, command = None):
        for ss in structures:
            spath = path + "/" + str(ss.sid)
            self.initdir(ss,spath)
            self.run_buildmatself(ss, spath)
        self.waitforproc(0.1)
        for ss in structures:
            spath = path + "/" + str(ss.sid)
            self.run_sewself(ss, spath, command)
            
    def run_buildmatself(self,structure,spath):
        su.dispatch(self.buildmatself, [self.structfilename], spath)
        
    def run_sewself(self,structure,spath,inputs = None):
        '''Run sewself for given structure and structure path (spath)
        inputs (Optional): List of sting of commands to be executed (a-w)
        '''
        commands = ""
        if inputs is None:
            inputs = self.commands
        
        for c in inputs:
            if c[0] == 'a':
                commands += "a\n" + c[1] + "\n" + c[2] + \
                 "\n" + c[3] + "\n"
            elif c[0] == 'b':
                commands += "b\n"
                commands += str(self.numpar.get('lattice_temp')) + "\n"
                commands += str(self.numpar.get('el_temp')) + "\n"
            elif c[0] == 'c':
                commands += "c\n"
            elif c[0] == 'd':
                commands += "d\n"
            elif c[0] == 'e':
                commands  += "e\n"   # compute self-consitent potential
                commands += str(self.numpar.get('lattice_temp')) + "\n"
                commands += str(self.numpar.get('el_temp')) + "\n"
            elif c[0] == 'f':
                commands  += "f\n"
            elif c[0] == 'g':
                commands  += "g\n"
                commands += str(self.numpar.get('el_temp')) + "\n"
            else:
                print("WARNING: Option " + c[0] + \
                " not implemented in Isewself!")
        commands += "q\n^C\n"
                
        dbg.debug("Running sewself:", callclass=self)
        dbg.flush()
        proc = su.dispatch(self.sewself, [], spath, 
                            infile = subprocess.PIPE,
                            outfile = subprocess.PIPE, 
                            errfile = subprocess.PIPE)
        dbg.debug("Communicating sewself: \n" + commands, callclass=self)
        dbg.flush()
        out, err = proc.communicate(commands)
        structure.output = out
            
    def initdir(self,structure,spath):
        su.mkdir(spath)
        self.writeParameterFile(spath)
        self.writeSewselfPar(spath)
        self.writeStructFile(structure, spath)
    
    def gatherResults(self, structures, path, pathresults = None, runprog=True):
        self.ebound = []
        self.dipoles = []
        with open(path+'/results.log','a') as f:
            f.write('# Results for structures:\nID | merit | \
            N times layer width | Ndop times (zi, zf, nvol) | N times x \n')
            for ss in structures:
                
                levels = self.readEbound(path + "/" + str(ss.dirname))
                dipoles = self.readDipoles(ss)
                ebounds = []
                [ebounds.append(float(l)) for l in levels]
                self.ebound.append(ebounds)
                ss.ebound = ebounds
                ss.dipoles = dipoles
                self.dipoles.append( self.readDipoles(ss) )
                
                f.write(str(ss.sid)+" ")
                f.write(str(self.getMerit(ss, path)) + " ")
                for layer in ss.layers:
                    f.write(str(layer.width)+" ")
                for doping in ss.dopings:
                    for val in doping:
                        f.write( str(val) +" ")
                for layer in ss.layers:
                    f.write(str(layer.material.x)+" ")
                
                              
                f.write("\n")
                
                if self.merit == self.merits['Chi2']:
                    with open(path+ ss.dirname+"/chi2.log", 'w') as chif:
                        domega = 0.001
                        gamma = self.target[2]
                        E1 = self.target[0]
                        chif.write('# E2, E3, |Chi2(E1,E2,E3=E1-E2)| ')
                        chif.write('with E1 fixed to E1 = ' + str(E1) + '\n')
                        chif.write("# gamma = " + str(gamma) + "\n")
                        
                        
                        for i in range(0, 1000):
                            E2 = i*domega
                            chif.write(str( E2 ) + " " + str( E1-E2 ) + " " )
                            chif.write(str( self.calcChi2(ss, E1, E2, gamma))
                                        + "\n")
                            
                
                
                
    def readEbound(self,path):
        with open(path + '/Ebound.dat','r') as f:
            line = f.read()
        return line.split()[1:]
    
    def readAbsorption(self,path):
        with open(path + '/absorption.dat') as f:
            line = f.read()
        return line.split('\n')[0:-1]
    
    def readDipoles(self, structure):
        out = structure.output.split()
        dipoles = []
        for i in range(len(out)):
            if out[i] == 'z':
                dipoles.append( float( out[i+2] ) )
        return dipoles
    
    def writeStructFile(self,structure,path):
        '''Writes the structure file self.structfilename.'''
        
        mat_list = []
        for layer in structure.layers:
            if layer.material not in mat_list:
                mat_list.append(layer.material)
        
        filepath = path + "/" + self.structfilename
        with open(filepath,'w') as f:
            f.write(str(self.numpar["efield"]) + " # efield\n")
            f.write(str(self.numpar["ldiff"]) + " # diffusion length\n")
            f.write(str(len(structure.layers)) + " # number of layers\n")
            f.write(str(self.numpar["nbswellinj"]) + " # nbswellinjector\n")
            f.write(str(self.numpar["nbswellact"]) + " # nbswellactive\n")
            f.write(str(self.numpar["coeffls"]) + " # coeff (longer/shorter structure)\n")
            f.write("switches for structure length\n")
            f.write(str(self.numpar["bool_inj"]) + " # injector only\n")
            f.write(str(self.numpar["bool_act"]) + " # active only\n")
            f.write(str(self.numpar["bool_one"]) + " # one period\n")
            f.write("material parameters\n")
            f.write(str(self.numpar["matchoice"]) + " matchoice\n")
            for mat in mat_list:
                f.write(str(mat) + " ")
                f.write(str(mat.params[mp.Ec]) + " ")
                if mat.x is None:
                    f.write(str(0))
                else:
                    f.write(str(mat.x))
                f.write("\n")
            f.write( "end discretisation of the potential\n" )
            f.write( str(self.numpar["inc0"]) + " inc0\n" )
            f.write( str(self.numpar["incw"]) + " incw\n" )
            f.write( "Structure (mat, thick (A), doping x 1e18 cm^-3) \n" )
            index = 0
            for layer in structure.layers:
                f.write(str(layer.material) + " " + str(10*layer.width) + " " )
                doping = 0
                l0 = structure.layerPos(index)
                l1 = structure.layerPos(index) + layer.width
                for dop in structure.dopings:
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
                f.write( str( doping*1e-18 ) + "\n")
                index += 1


    
    def writeParameterFile(self,path):
        '''Writes the parameter file "path/mat.par".'''
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
        '''Write the parameter file "path/sewself.par".'''
        
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
                
                
    def getMerit(self,structure,path):
        path = path + "/" + structure.dirname
        if self.merit == self.merits["absorption"]:
            line = self.readAbsorption(path)
            dE = 1e-3
            maxabs = 0
            for i in range(0,len(line)):
                f = float(line[i].split()[0])
                if f > float(self.target) - dE:
                    maxabs = max(float(line[i].split()[1]), maxabs)
                elif f > float(self.target) + dE:
                    break
            return maxabs
        
        elif self.merit == self.merits["QWIP"]:
            output = structure.output.split(' ')
            for il in range(0,len(output)):
                if output[il] == 'E21':
                    E21 = float(output[il + 2])
                    z = float(output[il + 6])
                    return -abs(E21 - float(self.target))*10 + z
            return 0
        
        elif self.merit == self.merits["Elase"]:
            output = structure.output.split(' ')
            for il in range(0,len(output)):
                if output[il] == 'Laser':
                    return abs(float(output[il+3]) - float(self.target))
            return 0
        elif self.merit == self.merits["DeltaE_12"]:
            levels = self.readEbound(path)
            ebound = []
            [ebound.append(float(l)) for l in levels]
            try:
                return -abs(ebound[1]-ebound[0] - float(self.target))
            except(ValueError,IndexError):
                return 'ERROR'
        elif self.merit == self.merits["Chi2"]:
            
            E1 = self.target[0]
            E2 = self.target[1]
            E3 = E1-E2
            gamma = self.target[2]
            
            
            elevel = structure.ebound
            dipoles = structure.dipoles
            dopdens = 0.
            for d in structure.dopings:
                dopdens += (d[1]-d[0])*d[2]
            dopdens /= structure.length
            
            if len(elevel) < 3 or len(dipoles) < 3:
                return 0
            
            chi2 = 1./(E1-elevel[2]+elevel[0] - gamma*1j)
            chi2 += 1./(-E2+elevel[1]-elevel[0] - gamma*1j)
            chi2 *= dipoles[0]*dipoles[1]*dipoles[2]/(E3-elevel[2]+elevel[1] - gamma*1j)
            chi2 *= -const.qe/const.eps0*dopdens*1e6*1e-30*1e12
            
            return np.abs( chi2 )
        else:
            print("Merit function " + str(self.merit) + " not implemented yet!")
            return "ERROR"
                
    def waitforproc(self,delay,message=None):
        pactive = True
        while pactive:
            if message is not None:
                print(message)
            pactive = False
            for p in self.processes:
                if self.pltfm.jobstatus(p):
                    pactive=True
                    #break
            time.sleep(delay)
            
    def calcChi2(self, structure, E1, E2, gamma):
        
        E3 = E1-E2
        
        elevel = structure.ebound
        dipoles = structure.dipoles
        dopdens = 0.
        for d in structure.dopings:
            dopdens += (d[1]-d[0])*d[2]
        dopdens /= structure.length
        
        if len(elevel) < 3 or len(dipoles) < 3:
            return 0
        
        chi2 = 1./(E1-elevel[2]+elevel[0] - gamma*1j)
        chi2 += 1./(-E2+elevel[1]-elevel[0] - gamma*1j)
        chi2 *= dipoles[0]*dipoles[1]*dipoles[2]/(E3-elevel[2]+elevel[1] - gamma*1j)
        chi2 *= -const.qe/const.eps0*dopdens*1e6*1e-30*1e12
        
        return np.abs( chi2 )
                        