'''
Created on 16 Mar 2018

@author: martin
'''

from interface import Interface
from structure.classes import MaterialPar as Par
from utils.systemutil import SystemUtil as su
from numerics.runplatf import Local
import time
from utils.debug import Debugger as dbg
import subprocess

class Inegf(Interface):
    '''
    classdocs
    '''
    
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
    
    # index of data in negft.dat
    idat = {"eFd":0,"omega":1,"eFacd":2,"j":3,"gain":4,"dk":5,"konv":6,"errdyn":7,"ierror":8}

    def __init__(self,binpath,pltfm,wellmaterial,einspath = "./"):
        '''
        Constructor
        '''
        super(Inegf,self).__init__(binpath, pltfm)
        self.numpar.update(self.negf_numpar)
        self.numpar["maxits"] = 40
        self.progwann = binpath+"wannier8.out"
        self.prognegft = binpath+"negft8mpi.out"
        self.proghdiag = binpath+"hdiag8.out"
        self.progbandplot = binpath+"bandplot8.out"
        self.wellmat = wellmaterial
        self.processes= []
        self.einspath = einspath
        
    def __str__(self):
        return "Inegf"
        
    def initdir(self,ss,path):
        pathNegf=path+"/IV/"
        su.mkdir(pathNegf)
        self.writeWannier(ss,path)
        self.writeMaterial(self.wellmat, "# test of function 29/1 2018",path)
        self.writeNegftInp(path, self.einspath ,pathNegf)
        
    def runStructures(self,structures,path):
        local = Local()
        for ss in structures:
            spath = path+"/"+str(ss.dirname)
            su.mkdir(spath)
            self.initdir(ss, spath)
            #proc = su.dispatch(self.progwann, "",spath)
            proc = self.pltfm.submitjob(self.progwann,[],spath,1,"00:10")
            self.processes.append(proc)
        dbg.debug("Starting Wannier program.....\n",dbg.verb_modes["verbose"],self)
        dbg.flush()
        # TODO do not wait for all processes-start each strucrure when ready!
        self.waitforproc(1)
        dbg.debug("Starting negf....\n",dbg.verb_modes["verbose"],self)
        dbg.flush()
        #del processes
        #processes = []
        for ss in structures:
            spath = path+"/"+str(ss.dirname)
            # replacing default value of Nper in scatt3.inp
            proc = local.submitandwait("sed",['-i',"--in-place=''",'2s/1/'+str(self.numpar["Nper"])+'/', "scatt3.inp"],spath)
            # to make sure file is closed:
            proc.communicate()
            proc = local.submitandwait("cp",["scatt3.inp","IV/"], spath)
            proc.communicate()
            
            proc = self.pltfm.submitjob(self.prognegft,[],spath+"/IV/")
            #proc = su.dispatch(self.prognegft, "",spath+"/IV/")
            self.processes.append(proc)
        return self.processes
        
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
        # close all finished processes:
        for p in self.processes:
            p.wait()
            del p
            
    def checkactive(self):
        pactive = False
        for p in self.processes:
                if self.pltfm.jobstatus(p):
                    pactive=True
        return pactive
        
    
    def gatherResults(self, structures, pathwd, pathresults = None):
        '''
        Write results to pathresults/results.log and run hdiag and bandplot
        in pathwd/i/IV/eins/x/ for each i and x. Stores WS resutls as a
        new attribute levels[directory][WS level][data field] in each 
        Structure object in the list structures.
        '''
        
        if(pathresults is None):
            pathresults = pathwd
        
        with open(pathresults+'/results.log','w') as f:
            f.write('# Results for structures:\nID | N times layer width | N times Mat | Merit\n')
            for ss in structures:
                ss.wslevels = []
                spath = pathwd + "/" + str(ss.dirname)
                try:
                    dirlist = su.listdirs(spath+"/IV/eins/")
                except (OSError, IOError):
                    print "WARNING: could not find directory: " + spath+"/IV/eins/"
                    continue
                dirs = dirlist[0].split()
                
                for folder in dirs:
                    einspath = spath+"/IV/eins/"+folder
                    omega0 = self.numpar["omega0"]
                    omegaf = self.numpar["domega"]
                    self.runHdiag(einspath,omega0=omega0,omegaf=omegaf, gamma=0.001)
                    self.runBandplot(einspath, ss)
                    ss.wslevels.append(self.getWSdata(einspath))
                
                f.write(str(ss.dirname)+" ")
                for layer in ss.layers:
                    f.write(str(layer.width)+" ")
                for layer in ss.layers:
                    f.write(str(layer.material.x)+" ")
                              
                f.write(str(self.getMerit(ss, pathwd)))
                f.write("\n")
                
                
    def runHdiag(self, path, zshift=0, omega0=0, omegaf=1, Nomega = 1000, Nper=1, gamma = 0.0001, Nk = 0, Ek = -1):
        if omega0 == -1:
            commands = str(zshift)+"\n-1\n0\n-1\n"
        else:
            commands = str(zshift)+"\n"+str(omega0)+"\n"+str(omegaf)+"\n"+str(Nomega)
            commands += "\n"+str(Nper) +"\n"+str(gamma)+"\n"+str(Nk)+"\n"+str(Ek)
        
        proc = su.dispatch(self.proghdiag, [], path,infile = subprocess.PIPE)
        proc.communicate(commands)
        su.waitforproc(proc, 0.1)
        return proc
    
    def getWSdata(self,path):
        levels = []
        with open(path+"/wslevels.dat") as ff:
            # read eFd
            line = ff.next().split()
            efd = line[2]
            # read 3 lines of comments
            for _ in range(0,3):
                ff.next()
            for _ in range(0,self.numpar["Nstates"]):
                line = ff.next().split()
                data = []
                data.append(path)
                data.append(float(efd))
                for element in line:
                    data.append(float(element))
                levels.append(data)
        return levels
        
    def runBandplot(self, path, structure, valence = 0, localdata = 1, eFd = 0, plotWS = 1, 
                    whichWS = 1, square = 1, renorm = 0.3, pmin = -1, pmax = 1, zmin = None, zmax = None):
        if zmin is None:
            zmin = -structure.length
        if zmax is None:
            zmax = 2*structure.length
        commands = str(valence) + "\n" + str(localdata) + "\n"
        if localdata == 1:
            commands += str(plotWS) + "\n"
            if plotWS == 1:
                commands+= str(whichWS) + "\n"
        elif localdata == 0:
            commands+= str(eFd) + "\n"
        else:
            raise Exception("Parameter not valid: plot local data = "+localdata)
        commands += str(zmin) + "\n" + str(zmax) + "\n" + str(pmin) + "\n" + str(pmax) + "\n" + str(square) + "\n" + str(renorm) + "\n"
        
        proc = su.dispatch(self.progbandplot, [], path, infile = subprocess.PIPE)
        proc.communicate(commands)
        return proc

    def getMerit(self,structure,path):
        path = path + "/" + structure.dirname + "/IV"
        
        results = []
        
        try:
            
            with open(path+"/negft.dat",'r') as f:
                for line in f:
                    if '#' in line or line.split()[Inegf.idat.get("ierror")] == '1'\
                            or line.split()[Inegf.idat.get("konv")] == 'NaN':
                        continue
                    results.append(line.split())
                if results == []:
                    return "NO CONV"
            f.closed
            
        except (OSError, IOError):
            print "\nWARNING: Could not find directory: " + path
            return "ERROR"
        except ValueError:
            print "\nWarning: Error when getting results from: " + path
            return "ERROR"
        
        values = []
        for i in range(0, len(Inegf.idat)):
            datval = []
            for row in results:
                datval.append(float(row[i]))
            values.append(datval)
        
        # get the gain of structure
        # TODO: now using last point, take maximum instead!
        if self.merit==Interface.merits.get("max gain") :

            out = max(values[Inegf.idat.get("gain")])
            
        elif self.merit == Interface.merits.get("(max gain)/(current density)"):
            
            gain = values[Inegf.idat.get("gain")]
            j = values[Inegf.idat.get("j")]
            gainj = []
            [gainj.append(gain[i]/j[i]) for i in range(0,len(gain))]
            out = max(gainj)
                
        elif self.merit == self.merits.get("estimated gain") :
            try:
                dirlist = su.listdirs(path+"/eins/")
            except (OSError, IOError):
                print "WARNING: could not find directory: " + path+"/IV/eins/"
                return "ERROR"
            dirs = dirlist[0].split()
            maxgain = []
            for folder in dirs:
                einspath = path+"/eins/"+folder
                with open(einspath+"/gainFGR.dat") as f:
                    for line in f:
                        linestr = line.split()
                        try:
                            maxgain.append( float( linestr[1] ) )
                        except( IndexError ):
                            pass
            out = max(maxgain)
            
                
        else:
            print "No such merit function!"
        
        return out
    
    def writeMaterial(self,material,nametag,dirpath = None):
        if dirpath is None:
            path="material.inp"
        else:
            path = dirpath + "/material.inp"
        with open(path, 'w') as f:
            f.write(nametag + "\n")
            f.write(str(material.params[Par.meff])+"\n")
            f.write(str(material.params[Par.eps0])+"\n")
            f.write(str(material.params[Par.epsinf])+"\n")
            f.write(str(material.params[Par.ELO])+"\n")
            f.write(str(material.params[Par.Vdef])+"\n")
            f.write(str(material.params[Par.vlong])+"\n")
            f.write(str(material.params[Par.massdens])+"\n")
            f.write(str(material.params[Par.molV]))
        f.closed

    def writeWannier(self,struct,dirpath=None):
        if dirpath is None:
            path="wannier8.inp"
        else:
            path = dirpath + "/wannier8.inp"
        
        # re-scaling of the CBO to lowest energy:
        Ec = []
        [Ec.append(l.material.params[Par.Ec]) for l in struct.layers]
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
                    f.write(str(l.material.params[Par.Ec]-Ecmin)+ " ")
                    f.write(str(l.material.params[Par.meff])+ " ")
                    f.write(str(l.material.params[Par.Valloy])+ " ")
                    f.write(" # Layer " + str(il))
                    il=il+1
                    f.write("\n")
                f.write(str(struct.layers[0].material.params[Par.Ep])+"\n") # TODO: take average Ep?
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
            print "WARNING: Directory "+dirpath+" not found!"
        
    def writeNegftInp(self,pathwann,patheins,dirpath=None):
        if self.numpar["boolEins"]:
            readeins = ".TRUE."
        else:
            readeins = ".FALSE."
        if dirpath is None:
            path = "negft7.inp"
        else:
            path = dirpath+"/negft7.inp"
        try:
            with open(path,"w") as f:
                f.write(str(self.numpar["Nper"])+" "+str(self.numpar["Nstates"])+" "+str(self.numpar["NE"]) + " " + str(self.numpar["Nk"]))
                f.write(" # Nper, Nnu, NE, Nk\n")
                f.write(str(self.numpar["Nh"])+" "+str(self.numpar["Igauge"]))
                f.write(" #\n")
                f.write(str(self.numpar["Tlattice"]))
                f.write(" #\n")
                f.write(str(self.numpar["efield0"])+" "+str(self.numpar["defield"])+" "+str(self.numpar["Nefield"])+" ")
                f.write(" #\n")
                f.write(str(self.numpar["efac0"])+" "+str(self.numpar["defac"])+" "+str(self.numpar["Nefac"])+" ")
                f.write(" #\n")
                f.write(str(self.numpar["omega0"])+" "+str(self.numpar["domega"])+" "+str(self.numpar["Nomega"])+" ")
                f.write(" #\n")
                f.write(str(self.numpar["Emaxadd"])+ " " + str(self.numpar["Eminsub"])+" #\n")
                f.write(str(self.numpar["gen"])+ " #\n")
                f.write(str(self.numpar["maxits"])+" "+str(self.numpar["Iconv"])+" "+str(self.numpar["Bei"])+" " \
                        + str(self.numpar["Nhist"]) + " #\n")
                f.write(pathwann + "/\n")
                f.write(readeins+"\n")
                f.write(patheins+"\n")
                if(self.numpar["use-poisson"]):
                    f.write(".TRUE.")
                else:
                    f.write(".FALSE.")
                f.write(" #\n")
                if(self.numpar["boolPrinc"]):
                    f.write(".TRUE.")
                else:
                    f.write(".FALSE.")
                f.write(" #\n")
                if(self.numpar["use-e-e"]):
                    f.write(".TRUE.")
                else:
                    f.write(".FALSE.")
                f.write(" #\n")
                
            f.closed
        except IOError:
            print "WARNING: Directory "+dirpath+" not found!"
    