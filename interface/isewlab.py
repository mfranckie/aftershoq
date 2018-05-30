'''
Created on 12 Mar 2018

@author: martin
'''

from interface import Interface
from structure.classes import MaterialPar as mp
from utils import const
from utils.systemutil import SystemUtil as su
import time

class Isewlab(Interface):
    '''
    Interface for sewself.
    '''
    # dictionary with name and default value
    
    numpar = {
        "efield0" : 0.0,
        "defield" : 0.0,
        "Nefield" : 0.0,
        "temp"   : 100,
        "ldiff"  : 0.0,
        "nbswellinj":4,
        "nbswellact":4,
        "coeffls":1,
        "bool_inj":0,
        "bool_act":0,
        "bool_one":1,
        "matchoice":1,
        "inc0":5,
        "incw":2,
        "verbosity" : "Verbose"
    }
    
    buildpot_params = {
        "bulk-step"               : 1,
        "interface-step"          : 0.01,
        "interface-diffusion"     : 0.0,  # DIFFUSION -- SEGREG
        "mesh-style"              : 'fixed-step',
        "phony-left-barrier"      : 'AUTO', # AUTO, USER, NO
        "phony-right-barrier"     : 'AUTO', #AUTO, USER, NO
        "has-box-wall"            : 'TRUE',# has box wall -- NEED TO BE TRUE
        "auto-box-wall"           : 'TRUE', # auto discontinuity
        
    }
        
    left_barrier = {
        
        "lb-thickness" : 350,
        "lb-material"  : 'AlGaAs',
        "lb-x" : 0.15,
        "lb-label" : '"left barrier"',
        "lb-mass" : 0.8,
        "lb-gap" : 0.1,
         "lb-discont" : 5.8,
    } 
        
    right_barrier = {
         
        "rb-thickness" : 35,
        "rb-material"  : 'AlGaAs',
        "rb-x" : 0.15,
        "rb-label" : '"right barrier"',
        "rb-mass" : 0.8,
        "rb-gap" : 0.1,
        "rb-discont" : 5.8,
    }
                       
    bw_layer = {
        
        "bw-thickness" : 2.0,
        "bw-material"  : 'AlAs',
        "bw-label" : '"Box Layer"'
    }
    
    solver_params = {
        
        "emin" : 0.001,
        "emax" : 1.6,
        "up-to-bound-state" : 20,
        "continuum-emin" : 0.001,
        "continuum-emax" : 2.0,
        "initial-samples"  : 1200,
        "np"  : 'TRUE',
        "energy-precision"        : 1e-14,
        "max-divergence"          : 1e-4,
        #//letal-divergence      = TRUE;
        "max-iterations"          : 10000,
        "boundify"                : 'TRUE',
        "adjust-box-wall"         : 'TRUE',
        "wf-resampling"           : 'TRUE',
        "wf-step"                 : 1.0,
        "up-to-continuum-state"  : 200
        
    }
    
    absorption_params = {
        
        "min-photon-energy" : 1e-5,
        "max-photon-energy" : 500e-3,
        "spectrum-sampling" : 512,
        "k-space-sampling" : 64,
        "number-of-kT-before-cut-off" : 7, # not too low
        "use-non-parabolicity" : 'TRUE',
        "default-subband-broadening" : 4.0e-3,  # 0.5*FWHM = upper-broadening + lower-broadening
        "modal-index" : 3.2,
        "modal-overlap" : 1.0      
        
    }
    
    ifr_params = {
        
        "angular-sampling" : 32,
        "ifr-inplane-corr" : 90,
        "ifr-vertical-corr" : 15.0,
        "ifr-height" : 1.2
        
    }
    
    hlo_params = {
        
       "level-of-details" : 8,
        "angular-sampling" : 32,
        "kp0" : 64,
        "exp-cutoff" : 1e-6
        
    }
    
    alloy_params = {
        
        "angular-sampling" : 32 # not used
        
    }

    imp_params = {
        
        "angular-sampling" : 32,
        "formfactor-lod" : 8,
        "exp-cutoff" : 1e-6,
        "crop-profile" : 'TRUE',
        "zero-trigger" : 0.0,
        "wf-step" : 3.0
        
    }
    
    transport_params = {
        
        "taup-sorting" : 0.08,
        "tunneltime-max" :  100,
        "tunneltime-maxstates" :  10,
        "initial-temperature" :  100.0, # Originally 1500 (necessary for the convergence)
        "k-space-sampling" :  256,
        "number-of-kT-before-cut-off" :  7,
        "hlo-energy" :  34.0e-3,
        "hlo-temperature" :  70.0, #300.0,
        "hlo-qscreen" :  0.0,
        "ifr-inplane-corr" :  90, #90
        "ifr-vertical-corr" :  15.0, #8, #30.0, #15.0,
        "ifr-height" :  1.2,#1.2
        "use-uniform-taup" :  'FALSE',
        "uniform-taup" :  0.04,
        "pop-tolerance" :  1e-3, # 1% Of ns
        "pop-max-iterations" :  500,
        "pop-stab-damping" :  1e-3,
        "electron-T-maxiter" :  1200,
        "electron-T-tolerance" :  0.5,     # in K
        "solution-maxiter" :  500,
        "current-uniformity-limit" :  1.0e-6,
        "solution-uniformity-limit" :  1.0e-2,
        "imaginary-part-limit" :  1.0e-10,

        # Light Parameters */
        "light-fixed-laser-energy" :  'FALSE',
        "light-use-bloch-gain" :  'TRUE',
        "light-laser-energy" :  50.0e-3,
        "light-gain-window-min-energy" :  40.0e-3,
        "light-gain-window-max-energy" :  60.0e-3,
        "light-gain-window-sampling"   :  64,
        "light-losses" :  72, # n655: 3mm HR (2.3 + 2.2+.5)/0.68 env 7cm-1
        "light-initial-photonflux" :  2e22,
        "light-bracketing-maxiter" :  128,
        "light-convergence-maxiter" :  5000,
        "light-convergence-tolerance" :  1e-1,     # in cm-1 */
        "light-damping-factor" :  0.3,
        "light-photonflux-precision" :  0.005
    }
    
    selfsolver_params = {
        
        'max-iterations' : 500,
        'damping-factor' : 1e-3,
        'energy-precision' : 1.5e-5,
        'period-wraping' : 2,
        'output-history-file' : 'FALSE',
        'convergence-crop' : 4
        
        }
    
    thermalmodel_params = {
        'initial-fermi-min'     : -1e-3,
        'initial-fermi-max'    : 1e-3,
        'fermi-bracketing-max-iterations' : 100,
        'fermi-brent-max-iterations' : 100,
        'fermi-tolerance' : 1e-8
    }
    
    show_options = {
        "underline-eigenstates" : 'TRUE'
    }



    def __init__(self, binpath, pltfm):
        '''
        Constructor
        '''
        super(Isewlab,self).__init__(binpath, pltfm)
        self.processes = []
        self.sewlab = binpath
        self.samplefilename = "structure.sample"
        self.structurename = 'AGS'
        self.scriptfilename = 'script.c'
        self.resultfile = "results.txt"
        self.popfile = 'pop.txt'
        self.dipolefile = 'dipoles.txt'
        self.ratefile = 'rates.txt'
        self.energiesfile = 'energies.txt'
        self.merits.update({'DeltaE_12' : 8, 'Elase' : 9})
        
    def __str__(self):
        return "Sewself"
        
    def runStructures(self,structures,path):
        for ss in structures:
            spath = path + "/" + str(ss.dirname)
            self.initdir(ss,spath)
            self.run_sewlab(ss, spath)
                
    def run_sewlab(self,structure,spath):
        proc = su.dispatch(self.sewlab, [self.scriptfilename], spath)
        self.processes.append(proc)
            
    def initdir(self,structure,spath):
        su.mkdir(spath)
        self.writeSampleFile(structure, spath)
        self.writeScriptFile(spath)
    
    def gatherResults(self,structures,path):
        for s in structures:
            s.results = self.readResults(s, path)
            s.populations = self.readPop(s, path)
            s.dipoles = self.readDipoles(s, path)
            s.energies = self.readEnergies(s, path)
            s.rates = self.readRates(s, path)
            
        with open(path+'/results.log','w') as f:
            f.write('# Results for structures:\nID | N times layer width | N times Mat | Merit\n')
            for ss in structures:
                f.write(str(ss.sid)+" ")
                for layer in ss.layers:
                    f.write(str(layer.width)+" ")
                for layer in ss.layers:
                    f.write(str(layer.material.x)+" ")
                              
                f.write(str(self.getMerit(ss, path)))
                f.write("\n")
            
    def readResults(self, s, path):
        results = []
        with open(path + s.dirname + "/" + self.resultfile, 'r') as f:
            for line in f:
                vals = []
                [vals.append(float(val)) for val in line.split()]
                results.append( vals )
        return results
    
    def readPop(self, s, path):
        results = []
        with open(path + s.dirname + "/" + self.popfile, 'r') as f:
            for line in f:
                vals = []
                [vals.append(float(val)) for val in line.split()]
                results.append( vals )
        return results
    
    def readDipoles(self, s, path):
        results = []
        with open(path + s.dirname + "/" + self.dipolefile, 'r') as f:
            for line in f:
                vals = []
                [vals.append(float(val)) for val in line.split()]
                results.append( vals )
        return results
    
    def readRates(self, s, path):
        results = []
        with open(path + s.dirname + "/" + self.ratefile, 'r') as f:
            for line in f:
                vals = []
                [vals.append(float(val)) for val in line.split()]
                results.append( vals )
        return results
    
    def readEnergies(self, s, path):
        results = []
        with open(path + s.dirname + "/" + self.energiesfile, 'r') as f:
            for line in f:
                results.append( float(line) )
        return results
    
    def writeSampleFile(self,structure,path):
        filepath = path + "/" + self.samplefilename
        with open(filepath,'w') as f:
            f.write('// ---- material system\n')
            f.write('materials {\n\n')
            materials = []
            matstr = []
            for l in structure.layers:
                if str(l.material) in matstr:
                    continue
                
                matstr.append(str(l.material))
                if l.material.x is None or l.material.x == 0 or l.material.x == 1:
                    self.writeBulk(l.material, f)
                    materials.append(l.material)
                else:
                    if str(l.material.mat1) not in matstr:
                        self.writeBulk(l.material.mat1,f)
                        matstr.append(str(l.material.mat1))
                        materials.append(l.material.mat1)
                    if str(l.material.mat2) not in matstr:
                        self.writeBulk(l.material.mat2,f)
                        matstr.append(str(l.material.mat2))
                        materials.append(l.material.mat2)
                    self.writeAlloy(l.material, f)

            for i in range(0, len(materials)-1 ):
                m1 = materials[i]
                for j in range (i+1,len(materials) ):
                    m2 = materials[j]
                    self.writeInterface(f, m2, m1)
                    
            f.write('}\n\n')
            
            f.write('// -----\n')
            f.write('// Definition of the MQW structure\n\n')
            f.write('sequence {\n')
            f.write('\tlabel = "'+self.structurename+'"; // Auto-Generated Structure\n')
            f.write('\txray = 1.0;\n')
            index = 0
            discont = 0
            for l in structure.layers:
                f.write('\tlayer { ')
                f.write('thickness = ' + str( l.width*10 ) + "; " )
                f.write('material = ' + str(l.material) + "; " )
                if l.material.x is None:
                    pass
                else:
                    f.write('x = ')
                    f.write(str(l.material.x))
                f.write("; ")
                f.write('mass = ' + str(l.material.params[mp.meff]) + "; " )
                f.write('gap = ' + str(l.material.params[mp.Eg]) + "; " )
                if index < len(structure.layers)-1:
                    discont = structure.layers[index + 1].material.params[mp.Ec]
                else:
                    discont = structure.layers[0].material.params[mp.Ec]
                discont = discont - l.material.params[mp.Ec]
                f.write('discont = ' + str(discont) + "; " )
                if abs(structure.layerPos(index)-structure.dopings[0][0]) < l.width:
                    f.write('doping = ' + str(structure.dopings[0][2]*1e-18) + ";")
                f.write(' }\n')
                index += 1
            f.write('}\n\n')
            
            f.write('// ----- buildpot numerical parameters\n')
            self.writeBuildPotParams(f)
            
            f.write('// ----- eigen solver parameters\n')
            f.write('solver {\n')
            for key,value in self.solver_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- absoption parameters\n')
            f.write('absorption { \n')
            for key,value in self.absorption_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- ifr parameters\n')
            f.write('ifr { \n')
            for key,value in self.ifr_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- hlo parameters\n')
            f.write('hlo { \n')
            for key,value in self.hlo_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
                
            f.write('// ----- alloy disorder parameters\n')
            f.write('alloy-disorder { \n')
            for key,value in self.alloy_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- impurities parameters\n')
            f.write('impurities { \n')
            for key,value in self.imp_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- transport parameters\n')
            f.write('transport { \n')
            for key,value in self.transport_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- self solver parameters\n')
            f.write('selfsolver { \n')
            for key,value in self.selfsolver_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- thermal model parameters\n')
            f.write('thermal-model { \n')
            for key,value in self.thermalmodel_params.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- show options\n')
            f.write('show-options { \n')
            for key,value in self.show_options.items():
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
    def writeBulk(self,material,f):
        f.write('\tbulk {\n')
        f.write('\t\talias\t= ')
        f.write(str(material) + ';\n')
        f.write('\t\tmass\t= ')
        f.write(str(material.params[mp.meff]) + ';\n')
        f.write('\t\tepsilon-inf\t= ')
        f.write(str(material.params[mp.epsinf]) + ';\n')
        f.write('\t\tepsilon-zero\t= ')
        f.write(str(material.params[mp.eps0]) + ';\n')
        f.write('\t\thLO\t= ')
        f.write(str(material.params[mp.ELO]) + ';\n')
        f.write('\t\thTO\t= ')
        # TODO: include ETO in material description
        f.write(str(material.params[mp.ELO]) + ';\n')
        f.write('\t\tdirect-gap\t= ')
        f.write(str(material.params[mp.Eg]) + ';\n')
        f.write('\t\tlattice-const\t= ')
        f.write(str(material.params[mp.lattconst]) + ';\n')
        f.write('\t}\n\n')
        
    def writeAlloy(self,material,f):
        f.write('\talloy {\n')
        f.write('\t\talias\t= ')
        f.write(str(material) + ';\n')
        f.write('\t\tzero-fraction = ')
        f.write(str(material.mat2) + ';\n')
        f.write('\t\tfull-fraction = ')
        f.write(str(material.mat1) + ';\n')
        f.write('\t}\n\n')
        
    def writeInterface(self,f,mat1,mat2):
        f.write('\tinterface {\n')
        f.write('\t\tleft-material = ')
        f.write(str(mat1) + ";\n")
        f.write('\t\tright-material = ')
        f.write(str(mat2) + ";\n")
        f.write('\t\tdiscontinuity = ')
        discont = mat2.params[mp.Ec]-mat1.params[mp.Ec]
        f.write(str(discont) + ";\n")
        f.write('\t}\n\n')
        
    def writeBuildPotParams(self,f):
        f.write("buildpot-params {\n\n")
        for key,value in self.buildpot_params.items():
            f.write('\t' + key + " = " + str(value)+";\n")
        f.write('\n\tleft-barrier {\n')
        for key, value in self.left_barrier.items():
            f.write('\t\t' + key[3:] + " = " + str(value)+";\n")
        f.write('\t}\n')
        f.write('\n\tright-barrier {\n')
        for key, value in self.right_barrier.items():
            f.write('\t\t' + key[3:] + " = " + str(value)+";\n")
        f.write('\t}\n')
        f.write('\n\tbox-wall-layer {\n')
        for key,value in self.bw_layer.items():
            f.write('\t\t' + key[3:] + " = " + str(value)+";\n")
        f.write('\t}\n')
        f.write('}\n')
            
            
    def writeScriptFile(self,path):
        filepath = path + "/" + self.scriptfilename
        with open(filepath,'w') as f:
            f.write('Verbosity ')
            f.write(self.numpar["verbosity"] + ';\n\n')
            f.write('// Structure And Parameters\n')
            f.write('mqw = (Load Sequence From "')
            f.write(self.samplefilename)
            f.write('" At "')
            f.write(self.structurename)
            f.write('");\n')
            f.write('params = (Load Tree From "')
            f.write(self.samplefilename)
            f.write('");\n\n')
            f.write('// Variables\n')
            f.write('efield = ')
            f.write(str(self.numpar['efield0']))
            f.write(';\n')
            f.write('temp = ')
            f.write(str(self.numpar['temp']))
            f.write(';\n')
            f.write('\n// Potential And Self Basis\n')
            f.write('pot = (Buildpot mqw Using params);\n')
            f.write( "Sweep efield from ")
            f.write( str( self.numpar['efield0'] ) + " to " )
            f.write( str( self.numpar['Nefield']*self.numpar['defield']) )
            f.write(" step " + str( self.numpar["defield"] ) + "\n" )
            f.write('bpot = (Bias pot To efield);\n')
            f.write('sol = (Selftransport bpot Using params --no-superself --no-kinetic-balance)\n')
            f.write('Save sol.basis "wavef.txt":"wf"\n')
            f.write('sol.Lifetimes\n')
            f.write('Write efield sol.J (Stats Max sol.GainLorentz) (Stats MaxLoc sol.GainLorentz) ')
            f.write('To "')
            f.write(self.resultfile)
            f.write('"\n')
            f.write('Save sol "pop.txt":"populations"\n')
            f.write('Save sol "rates.txt":"rates"\n')
            f.write('Save sol "energies.txt":"energies"\n')
            f.write('Save sol "dipoles.txt":"dipoles"\n')
            f.write('Save sol.GainLorentz "gain.itx"\n')
            f.write("Endsweep\n")
            f.write("Quit\n")
                
    def getMerit(self,structure,path):
        path = path + "/" + structure.dirname
        if self.merit == self.merits["Elase"]:
            proc = su.dispatch('tail', ['-n','1',self.resultfile], path)
            [output, _] = proc.communicate()
            output = output.split()
            return abs(float(output[2]) - float(self.target))
        elif self.merit == self.merits["max gain"]:
            gain = []
            [gain.append(float(l[2])) for l in structure.results]
            return max(gain)
        else:
            print "Merit function " + str(self.merit) + "not implemented yet!"
            return "ERROR"
                
    def waitforproc(self,delay,message=None):
        pactive = True
        while pactive:
            if message is not None:
                print message
            pactive = False
            for p in self.processes:
                if self.pltfm.jobstatus(p):
                    pactive=True
                    #break
            time.sleep(delay)
                        