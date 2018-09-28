'''
Created on 12 Mar 2018

@author: Martin Franckie
'''

from .interface import Interface
import structure.matpar as mp
from utils import const
import utils.systemutil as su
import time
import numpy as np
import utils.debug as dbg

class Isewlab(Interface):
    '''Interface for sewlab, which is documented in:
    [Terazzi2010] Terazzi and Faist, N. J. Phys. 12, 033045 (2010)
    [Terazzi2011] Terazzi, PhD thesis, ETH Zuerich (2011)
    '''
    
    # dictionary with name and default value
    
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
        "light-losses" :  72,
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
    
    selftransport_flags = {
        "no-kinetic-balance" : True,
        "no-superself" : True,
        "compute-light" : False,
        "direct-scattering" : False,
        "edge-gain-line" : False,
        "no-ando-dephasing" : False
    }
    
    script_params = {
        "split-pot-layer": 1,
        "split-pot-bestFGR": True
    }
    
    light_params = {
        "width" : 4, # ridge width in micron
        "mirror_loss" : 7.7, # cm^-1
        "periods" : 40, # number of periods
        "Rfacet" : 0.25 # facet reflectivity
    }



    def __init__(self, binpath, pltfm, wellmaterial, version = '4.6.4'):
        '''
        Constructor. Subcalss-specific parameters:
        wellmaterial: Material of the well of the structure (used for
        dielectric and phonon properties).
        vesrion: Defines the version of the sewlab binary which is to be used.
        '''
        
        super(Isewlab,self).__init__(binpath, pltfm)
        self.numpar.update({ 
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
            })
        self.version = version
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
        self.wavefile = 'wavef.txt'
        self.potfile = 'pot.itx'
        self.dopingfile = 'doping.itx'
        self.bandplotfile = 'bandplot.txt'
        self.merits.update({
            'DeltaE_12' : 8, 'Elase' : 9
            })
        self.transport_params["hlo-energy"] = wellmaterial.params[mp.ELO]
        
    def __str__(self):
        return "Sewself"
        
    def runStructures(self,structures,path):
        '''Run simulations for all structures in the given structure list with
        the base path "path". This method dispatches all processes and returns
        the user has to wait for processes to finish before accessing results.
        
        Stores started processes in self.processes
        '''
        
        for ss in structures:
            spath = path + "/" + str(ss.dirname)
            su.mkdir(spath)
            spath = spath + "/" + str( abs( self.numpar["efield0"] ) ) + "/" 
            self.initdir(ss,spath)
            self.run_sewlab(ss, spath)
                
    def run_sewlab(self,structure,spath):
        
        proc = self.pltfm.submitjob(self.sewlab,[self.scriptfilename],spath)
        
        self.processes.append(proc)
            
    def initdir(self,structure,spath):
        '''Initialize the dierctory for Structure s, with base path "path".'''
        su.mkdir(spath)
        self.writeSampleFile(structure, spath)
        self.writeScriptFile(spath)
    
    def gatherResults(self,structures,path, wavescale = 1, 
                      square = True, pathresults=None, runprog=None):
        '''Write results to pathresults/results.log. Stores results in local
        variables results, populations, dipoles, energies, and rates. Also
        writes the wave functions to disc.
        '''
        
        for s in structures:
            s.results = self.readResults(s, path)
            s.populations = self.readPop(s, path)
            s.dipoles = self.readDipoles(s, path)
            s.energies = self.readEnergies(s, path)
            s.rates = self.readRates(s, path)
            self.saveBands(s, path, wavescale, square)
            
        with open(path+'/results.log','a') as f:
            f.write('# Results for structures:\nID | Merit | N times layer width | N times Mat \n')
            for ss in structures:
                f.write(str(ss.sid)+" ")
                f.write(str(self.getMerit(ss, path))+" ")
                for layer in ss.layers:
                    f.write(str(layer.width)+" ")
                for layer in ss.layers:
                    f.write(str(layer.material.x)+" ")
                              
                
                f.write("\n")
            
    def saveBands(self, s, path, wavescale, square):
        '''Read wave functions and potential profile from saved files and
        save them in self.bandplotfile for later plotting.
        '''
        
        spath = path + "/" + s.dirname
        if self.version == '4.6.4':
            splitchar = ","
        elif self.version == '4.6.5':
            splitchar = " "
        else:
            splitchar = ","
        for dir in su.listdirs(spath):
            try:
                fwave = open(spath + "/" + dir + "/" + self.wavefile, 'r')
                fpot = open(spath + "/" + dir + "/" + self.potfile, 'r')
                fbands = open(spath + "/" + dir + "/" + self.bandplotfile, 'w')
                
                for _ in range(0, 3):
                    fpot.readline()
                E = fwave.readline().split(splitchar)
                
                for lwave in fwave:
                    lpot = fpot.readline().split()
                    lwavesplit = lwave.split(splitchar)
                    fbands.write(lwavesplit[0] + " " + lpot[1])
                    for iwave in range(1, len(E)):
                        if square:
                            value = float( lwavesplit[iwave] )*float( lwavesplit[iwave] )*wavescale \
                                + float(E[iwave])
                        else:
                            value = float( lwavesplit[iwave] )*wavescale + float(E[iwave])
                        fbands.write(" " + str( value ))
                    fbands.write("\n")
                    
                fwave.close()
                fpot.close()
                fbands.close()
            except(IOError):
                pass
            
            
    def readResults(self, s, path):
        '''Read and return results for Structure s with base path "path".
        The results are read from self.resultfile and the returned results
        are sorted according to increasing electric field.
        '''
        
        results = []
        spath = path + "/" + s.dirname
        for dir in su.listdirs(spath):
            try:
                with open(spath + "/" + dir + "/" + self.resultfile, 'r') as f:
                    for line in f:
                        vals = []
                        [vals.append(float(val)) for val in line.split()]
                        results.append( vals )
            except(IOError):
                dbg.debug("WARNING: Error in directory: " + spath + "/" + dir, 
                          callclass = self)
        
        if results == []:
            return "ERROR"        
        # sort the results according to efield:
        results = np.array(results)
        k = results[results[:,0].argsort()]
        return k
    
    def readPop(self, s, path):
        '''Read and return populations for all levels for Structure s with
        base path "path".
        '''
        
        results = []
        spath = path + "/" + s.dirname
        for dir in su.listdirs(spath):
            tmp = []
            try:
                with open(spath + "/" + dir + "/" + self.popfile, 'r') as f:
                    for line in f:
                        vals = []
                        [vals.append(float(val)) for val in line.split()]
                        tmp.append( vals )
                results.append( tmp )
            except(IOError):
                dbg.debug("WARNING: Error in directory: " + spath + "/" + dir,
                          callclass = self)
        return results
    
    def readDipoles(self, s, path):
        '''Read and return dipoles between all levels for Structure s with
        base path "path".
        '''
        
        results = []
        spath = path + "/" + s.dirname
        for dir in su.listdirs(spath):
            tmp = []
            try:
                with open(spath + "/" + dir + "/" + self.dipolefile, 'r') as f:
                    for line in f:
                        vals = []
                        [vals.append(float(val)) for val in line.split()]
                        tmp.append( vals )
                results.append( tmp )
            except(IOError):
                print("WARNING: Error in directory: " + spath + "/" + dir)
        return results
    
    def readRates(self, s, path):
        '''Read and return scattering rates for all levels for Structure s
        with base path "path".
        '''
        
        results = []
        spath = path + "/" + s.dirname
        for dir in su.listdirs(spath):
            tmp = []
            try:
                with open(spath + "/" + dir + "/" + self.ratefile, 'r') as f:
                    for line in f:
                        vals = []
                        [vals.append(float(val)) for val in line.split()]
                        tmp.append( vals )
                results.append( tmp )
            except(IOError):
                print("WARNING: Error in directory: " + spath + "/" + dir)
        return results
    
    def readEnergies(self, s, path):
        '''Read and return energies for all levels for Structure s with
        base path "path".
        '''
        
        results = []
        spath = path + "/" + s.dirname
        for dir in su.listdirs(spath):
            tmp = []
            try:
                with open(spath + "/" + dir + "/" + self.energiesfile, 'r') as f:
                    for line in f:
                        tmp.append( float(line) )
                results.append( tmp )
            except(IOError):
                print("WARNING: Error in directory: " + spath + "/" + dir)
        return results
    
    def writeSampleFile(self,structure,path):
        '''Write the sample input file.'''
        
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
                
                doping = structure.layerDoping2D(index)*1e-18
                if doping > 0.0:
                    f.write('doping = ' + str(doping) + ";")
                
                f.write(' }\n')
                index += 1
            f.write('}\n\n')
            
            f.write('// ----- buildpot numerical parameters\n')
            self.writeBuildPotParams(f)
            
            f.write('// ----- eigen solver parameters\n')
            f.write('solver {\n')
            for key,value in list(self.solver_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- absoption parameters\n')
            f.write('absorption { \n')
            for key,value in list(self.absorption_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- ifr parameters\n')
            f.write('ifr { \n')
            for key,value in list(self.ifr_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- hlo parameters\n')
            f.write('hlo { \n')
            for key,value in list(self.hlo_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
                
            f.write('// ----- alloy disorder parameters\n')
            f.write('alloy-disorder { \n')
            for key,value in list(self.alloy_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- impurities parameters\n')
            f.write('impurities { \n')
            for key,value in list(self.imp_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- transport parameters\n')
            f.write('transport { \n')
            for key,value in list(self.transport_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- self solver parameters\n')
            f.write('selfsolver { \n')
            for key,value in list(self.selfsolver_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- thermal model parameters\n')
            f.write('thermal-model { \n')
            for key,value in list(self.thermalmodel_params.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
            f.write('// ----- show options\n')
            f.write('show-options { \n')
            for key,value in list(self.show_options.items()):
                f.write('\t' + key + " = " + str(value) + ";\n")
            f.write('}\n\n')
            
    def writeBulk(self,material,f):
        '''Write a bulk Material "material" to file "f".'''
        
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
        '''Write an alloy Material "material" to file "f".'''
        
        f.write('\talloy {\n')
        f.write('\t\talias\t= ')
        f.write(str(material) + ';\n')
        f.write('\t\tzero-fraction = ')
        f.write(str(material.mat2) + ';\n')
        f.write('\t\tfull-fraction = ')
        f.write(str(material.mat1) + ';\n')
        f.write('\t}\n\n')
        
    def writeInterface(self,f,mat1,mat2):
        '''Write an interface between Materials "mat1" and "mat2" 
        to file "f".
        '''
        
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
        '''Write buildot params to file "f".'''
        
        f.write("buildpot-params {\n\n")
        for key,value in list(self.buildpot_params.items()):
            f.write('\t' + key + " = " + str(value)+";\n")
        f.write('\n\tleft-barrier {\n')
        for key, value in list(self.left_barrier.items()):
            f.write('\t\t' + key[3:] + " = " + str(value)+";\n")
        f.write('\t}\n')
        f.write('\n\tright-barrier {\n')
        for key, value in list(self.right_barrier.items()):
            f.write('\t\t' + key[3:] + " = " + str(value)+";\n")
        f.write('\t}\n')
        f.write('\n\tbox-wall-layer {\n')
        for key,value in list(self.bw_layer.items()):
            f.write('\t\t' + key[3:] + " = " + str(value)+";\n")
        f.write('\t}\n')
        f.write('}\n')
            
            
    def writeScriptFile(self,path):
        '''Write the script file executed by sewlab from working directory
        "path".'''
        
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
            f.write(str(self.numpar['Tlattice']))
            f.write(';\n')
            f.write('\n// Potential And Self Basis\n')
            f.write('pot = (Buildpot mqw Using params);\n')
            f.write('d = (Period pot);\n')
            # If you want to sweep
            if( self.numpar['Nefield'] > 1 ):
                f.write( "Sweep efield from ")
                f.write( str( self.numpar['efield0'] ) + " to " )
                f.write( str( self.numpar['efield0'] + 
                              (self.numpar['Nefield']-1)*self.numpar['defield']) )
                f.write(" step " + str( self.numpar["defield"] ) + "\n" )
            f.write('bpot = (Bias pot To (* efield d) );\n')
            
            if(self.script_params["split-pot-layer"] > -1):
                if(self.script_params["split-pot-bestFGR"] == True):
                    f.write('layer = (BestFGR bpot 1 Using params);\n')
                else:
                    f.write('layer = ' + 
                            str( self.script_params["split-pot-layer"]) + ';\n' )
                f.write('left = (Splitpot bpot At-layer layer And-get Left);\n')
                f.write('right = (Splitpot bpot At-layer layer And-get Right);\n')
                f.write('sol = (Selftransport left right Using params')
        
            else:
                f.write('sol = (Selftransport bpot Using params')
                
            for flag in list(self.selftransport_flags.items()):
                if flag[1] is True:
                    f.write(' --' + flag[0])
            f.write(')\n')
            f.write('Save sol.basis "' + self.wavefile + '":"wf"\n')
            f.write('Save bpot "' + self.potfile + '"\n')
            f.write('Write efield sol.J (Stats Max sol.GainLorentz) (Stats MaxLoc sol.GainLorentz) ')
            if self.selftransport_flags["compute-light"] == True:
                f.write(' (OpticalPower sol ')
                f.write(str( self.light_params["periods"] ) + ' ')
                f.write(str( self.light_params["width"] ) +  ' ')
                f.write(str( self.light_params["Rfacet"] ))
                f.write(' ) (Wallplug sol ')
                f.write(str( self.light_params["mirror_loss"] ))
                f.write(') ')
            f.write('To "')
            f.write(self.resultfile)
            f.write('"\n')
            f.write('Save sol "' + self.popfile + '":"populations"\n')
            f.write('Save sol "' + self.ratefile + '":"rates"\n')
            f.write('Save sol "' + self.energiesfile + '":"energies"\n')
            f.write('Save sol "' + self.dipolefile + '":"dipoles"\n')
            f.write('Save bpot.doping "' + self.dopingfile + '"\n')
            f.write('Save sol.GainLorentz "gain"!efield!".itx"\n')
            if( self.numpar['Nefield'] > 1 ):
                f.write("Endsweep\n")
            f.write("Quit\n")
                
    def getMerit(self,structure,path):
        '''Returns the merit function evaluated for the Structure structure,
        with base path "path". 
        '''
        
        path = path + "/" + structure.dirname
        try:
            
            if self.merit == self.merits["Elase"]:
                proc = su.dispatch('tail', ['-n','1',self.resultfile], path)
                [output, _] = proc.communicate()
                output = output.split()
                return abs(float(output[2]) - float(self.target))
            elif self.merit == self.merits["max gain"]:
                gain = []
                [gain.append(float(l[2])) for l in structure.results]
                return max(gain)
            elif self.merit == self.merits["wall plug efficiency"]:
                wp = []
                if structure.results == "ERROR":
                    return "ERROR"
                [wp.append(float(l[5])) for l in structure.results]
                return max(wp)
            else:
                print("Merit function " + str(self.merit) + "not implemented yet!")
                return "ERROR"
        
        # Catch index errors arising from sewlab failing
        except(IndexError):
            return "ERROR"
                
            
    def setTe(self, Te):
        self.Te = Te
        self.numpar["Te"] = Te
        self.transport_params["initial-temperature"] = Te
        
    def setTlattice(self, TL):
        self.TL = TL
        self.transport_params["hlo-temperature"] = TL
        
    def useKinBal(self, bool):
        if( bool == True ):
            self.script_params["no-kinetic-balance"] = False
        else:
            self.script_params["no-kinetic-balance"] = True
            
    def useSuperself(self, bool):
        if( bool == True ):
            self.script_params["no-superself"] = False
        else:
            self.script_params["no-superself"] = True
            
    def splitPot(self, default_layer, use_bestFGR = True):
        '''Use BestFGR to split potentail at default_layer. If
        use_bestFGR = True (default) then sewlab uses BestFGR to split 
        the potential at the best layer for FGR calculations in the two
        regions. The splitting is done for each sweep point.
        Splitting can be turned off by setting default_layer = -1.
        '''
        
        self.script_params["split-pot-layer"] = default_layer
        self.script_params["split-pot-bestFGR"] = use_bestFGR
        
        
                        
    def computeLight(self, bool, minE = None, maxE = None, wg_loss = 4., m_loss = 4., width = 4, Rfacet = 0.25, periods = 30):
        '''Choose to compute light. Optional parameters if bool = True:
        minE/maxE: min/max photon energy to evaluate photon flux
        wg_loss/m_loss: total waveguide / mirror losses
        width: Ridge width in micrometers
        Rfacet: Facet reflectivity (of one facet)
        periods: Number of periods of the active region
        '''
        if( bool == True ):
            self.script_params["compute-light"] = True
            self.transport_params["light-losses"] = wg_loss
            self.transport_params["light-gain-window-min-energy"] = minE
            self.transport_params["light-gain-window-max-energy"] = maxE
            self.light_params["mirror_loss"] = m_loss
            self.light_params["Rfacet"] = Rfacet
            self.light_params["width"] = width
            self.light_params["periods"] = periods
        else:
            self.script_params["compute-light"] = False