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
from scipy.interpolate import interp1d
import threading
import os.path

class Inegf(Interface):
    '''
    Inteface for the NEGF8 program package, documented in:
    [Wacker2013] Wacker et al., IEEE Sel. Top. Quant. Electron. 19,
                 1200611 (2013)
    [Lindskog2013] Lindskog et al., Proc. SPIE 8846, 884603 (2013)
    [Franckie2016] Franckie, PhD Thesis, Lund Univesity (May, 2016)
    [Winge2016a]   Winge et al., J. Phys.: Conf. Series 696,
                 012013 (2016)
    [Winge2016b]   Winge, PhD Thesis, Lund University (Nov, 2016)
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

    hdiag_numpar = {
        "fgr_omega0" : 0,
        "fgr_omegaf" : 0.1,
        "fgr_Nomega" : 1000,
        "fgr_gamma"  : 0.005,
        "zshift_trials" : 2
        }

    # index of data in negft.dat
    idat = {"eFd":0,"omega":1,"eFacd":2,"j":3,"gain":4,"dk":5,"konv":6,"errdyn":7,"ierror":8}

    def __init__(self,binpath="",pltfm=Local(),wellmaterial=GaAs(),einspath = "./"):
        '''
        Constructor. Subclass specific parameters:
        wellmaterial: Material of the well of the structure (used for
        dielectric and phonon properties).
        einspath: Path to eins-files used for accelerating convergence.
        '''

        super(Inegf,self).__init__(binpath, pltfm)
        self.numpar.update(self.negf_numpar)
        self.numpar.update(self.hdiag_numpar)
        self.numpar["maxits"] = 40
        self.progwann = binpath+"wannier8.out"

        if pltfm.paral == pltfm.paral_modes["SERIAL"]:
            self.prognegft = binpath + "negft8.out"
        elif pltfm.paral == pltfm.paral_modes["OMP"]:
            self.prognegft = binpath + "negft8omp.out"
        elif pltfm.paral == pltfm.paral_modes["MPI"]:
            self.prognegft = binpath + "negft8mpi.out"

        self.proghdiag = binpath+"hdiag8.out"
        self.progbandplot = binpath+"bandplot8.out"
        self.wellmat = wellmaterial
        self.processes= []
        self.einspath = einspath
        self.datpath = "Run"
        self.target = None

    def __str__(self):
        return "Inegf"

    def initdir(self,ss,path):
        '''Initialize the dierctory for Structure s, with base path "path".'''
        pathNegf=os.path.join(path,self.datpath)
        su.mkdir(pathNegf)
        self.writeWannier(ss,path)
        self.writeMaterial(self.wellmat, "# "+str(self.wellmat.name),path)
        self.writeNegftInp(su.abspath(path), self.einspath ,pathNegf)

    def runStructures(self, structures, path, runwannier = True):
        '''Run simulations for all structures in the given structure list with
        the base path "path". This method dispatches all processes and returns
        the user has to wait for processes to finish before accessing results.

        Stores started processes in self.processes
        '''

        local = Local()
        for ss in structures:
            spath = path+"/"+str(ss.dirname)
            su.mkdir(spath)
            self.initdir(ss, spath)
            if runwannier:
                proc = self.pltfm.submitjob(self.progwann,[],spath,1,"00:10")
                self.processes.append(proc)
        if runwannier:
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
            if not runwannier:
                self.writeNegftInp(su.abspath(spath), self.einspath ,spath+"/"+self.datpath)
            # replacing default value of Nper in scatt3.inp
            proc = local.submitandwait("sed",['-i',"--in-place=''",'2s/1/'+str(self.numpar["Nper"])+'/', "scatt3.inp"],spath)
            # to make sure file is closed:
            out, err = proc.communicate()
            proc = local.submitandwait("cp",["scatt3.inp","gw.inp",os.path.join(".",self.datpath)], spath)
            out, err = proc.communicate()

            proc = self.pltfm.submitjob(self.prognegft,[],os.path.join(spath, self.datpath))

            self.processes.append(proc)
        return self.processes


    def runStructSeq(self, structures, path, seq = None, runwannier = True):
        """Run NEGF sequentially, two times with (possibly) different parameters.
        Generates one Thread for each structure.

        Parameters:

        structures: list[Structure]
            Structures to evaluate. The same evaluation will be carried out for
            each Structure object in the list.
        path: String
            Base path to directory with structure sub-directories.
        seq: list[dictionary]
            List with two dictionaries of the form of Inegf.numpar, that will
            be used in order.
        runwannier : boolean
            Specifies whether the wannier program will be run. Defaults to True.

        Returns: list[Thread]
            A list with the started threads. The calling function may implement

                for tt in threads: tt.join()

            in order to wait for all threads to complete.
        """

        threads = []
        for ss in structures:
            t = threading.Thread(target=self.runSequence
                                 , args = (ss, path, seq, runwannier) )
            t.start()
            threads.append(t)

        return threads

    def runSequence(self, structure, path, seq = None, runwannier = True):
        """Run NEGF sequentially for a single structure on a single Thread.

        Parameters:

        structures: list[Structure]
            Structures to evaluate. The same evaluation will be carried out for
            each Structure object in the list.
        path: String
            Base path to directory with structure sub-directories.
        seq: list[dictionary]
            List with two dictionaries of the form of Inegf.numpar, that will
            be used in order.
        runwannier : boolean
            Specifies whether the wannier program will be run. Defaults to True.

        Returns: numpy.Array, numpy.Array
            negft_iv, negft_gain contains the results for the IV and gain simulations
            on matrix form.

        """
        if seq is None:
            numpar = self.numpar.copy
            seq = [numpar,numpar]

        spath = os.path.join(path,str(structure.dirname))
        su.mkdir(spath)
        self.initdir(structure, spath)
        if runwannier:
            # Wannier program
            proc = self.pltfm.submitjob(self.progwann, [], spath, 1, "00:10")
            print(f"sid={structure.sid} running Wannier...")
            while self.pltfm.jobstatus(proc):
                time.sleep(1)

        self.runNEGF(spath, numpar=seq[0], datpath= "IV")
        
        try:
            negft_iv = self.getresults(structure, path, datpath="IV")
        except FileNotFoundError as e:
            print("Simulation failed: ")
            print(e)
            return False

        # interpolate to find maximum IV point
        x = np.linspace(negft_iv[0,0],negft_iv[-1,0])
        f = interp1d(negft_iv[:,0], negft_iv[:,3], kind='cubic')
        imax = np.argmax(f(x))
        # closest point for which we have ein files
        i = np.argmax(negft_iv[:,3])
        print(f"sid={structure.sid}: jmax = {negft_iv[i,3]}, vmax = {negft_iv[i,0]}")

        numpar = seq[1]

        numpar["efield0"] = x[imax]

        einspath = f"{negft_iv[i,0]*1000:.1f}_{negft_iv[i,1]*1000:.1f}_{negft_iv[i,2]*1000:.1f}/"
        einspath=os.path.join("..","IV","eins",einspath)

        self.runNEGF(spath, numpar=numpar, einspath=einspath)

        negft_gain = self.getresults(structure, path)

        return negft_iv, negft_gain

    def runNEGF(self, spath, einspath=None, datpath=None, numpar=None):
        """Start the NEGF program in a specified directory.

        Parameters:
        spath: String
            Path to structure directory where execution should start.
        einspath: String
            (Optional) Relative path to eins files.
        datpath: String
            (Optional) Name of the execution directory. Defaults to self.datpath
        numpar: dictionary
            (Optional) Dictionary with numerical parameters. Defaults to self.numpar.
        """

        if einspath is None: einspath = self.einspath
        if datpath is None: datpath = self.datpath
        if numpar is None: numpar = self.numpar

        local = Local()
        su.mkdir(os.path.join(spath,datpath))

        self.writeNegftInp(su.abspath(spath), einspath ,
                           os.path.join(spath,datpath), numpar=numpar)
        # replacing default value of Nper in scatt3.inp
        proc = local.submitandwait("sed",['-i',"--in-place=''",
                                          '2s/1/'+str(numpar["Nper"])+'/',
                                          "scatt3.inp"],spath)
        # to make sure file is closed:
        out, err = proc.communicate()
        proc = local.submitandwait("cp",["scatt3.inp","gw.inp",
                                         os.path.join(".",datpath)], spath)
        out, err = proc.communicate()

        proc = self.pltfm.submitjob(self.prognegft,[],os.path.join(spath,datpath))

        print(f"Running NEGF in {spath}/{datpath}...")
        while self.pltfm.jobstatus(proc):
            time.sleep(1)
        print(f"Finished NEGF in {spath}/{datpath}!")

    def checkactive(self):
        pactive = False
        for p in self.processes:
                if self.pltfm.jobstatus(p):
                    pactive=True
        return pactive


    def gatherResults(self, structures, pathwd, pathresults = None, runhdiag = True, runbandplot = True):
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
                ss.wslevels = []
                ss.dipoles = []
                spath = pathwd + "/" + str(ss.dirname)
                try:
                    dirlist = su.listdirs(os.path.join(spath,self.datpath,"eins"))
                except (OSError, IOError):
                    print("WARNING: could not find directory: " + os.path.join(spath,self.datpath,"eins"))
                    continue
                dirs = dirlist
                for folder in dirs:
                    einspath = os.path.join(spath,self.datpath,"eins",folder)
                    omega0 = self.numpar["fgr_omega0"]
                    omegaf = self.numpar["fgr_omegaf"]
                    Nomega = self.numpar["fgr_Nomega"]
                    gamma  = self.numpar["fgr_gamma"]
                    if runhdiag:
                        self.runHdiag(einspath,omega0=omega0,omegaf=omegaf, Nomega = Nomega, gamma=gamma)

                        # Check if daigolalization failed, try zshift_trial times with zshift:
                        zshift = 0.
                        for _ in range(self.hdiag_numpar["zshift_trials"]-1):
                            if self.checkWSdens(einspath) == False:
                                zshift += ss.length/self.hdiag_numpar["zshift_trials"]
                                self.runHdiag(einspath,zshift=zshift, omega0=omega0,omegaf=omegaf,
                                          Nomega = Nomega, gamma=gamma)
                            else:
                                break
                        if runbandplot:
                            self.runBandplot(einspath, ss)

                    (levels, dipoles) = self.getWSdata(einspath)
                    ss.wslevels.append(levels)
                    ss.dipoles.append(dipoles)

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
                            domega = 0.001
                            if len( self.target ) > 0:
                                E1 = self.target[0]
                            else:
                                E1  = ss.wslevels[0][2][4]-ss.wslevels[0][0][4]
                            if len( self.target ) < 3:
                                gamma = None
                            else:
                                gamma = self.target[2]
                            chif.write('# E2, E3, |Chi2(E1,E2,E3=E1-E2)| ')
                            chif.write('with E1 fixed to E1 = ' + str(E1) + '\n')
                            if gamma is None:
                                chif.write("# gamma from NEGF \n")
                            else:
                                chif.write("# gamma = " + str(gamma) + "\n")

                            for i in range(0, 1000):
                                E2 = i*domega
                                chif.write(str( E2 ) + " " + str( E1-E2 ) + " " )
                                chif.write(str( self.calcChi2(ss, E1, E2, gamma))
                                            + "\n")


    def runHdiag(self, path, zshift=0, omega0=0, omegaf=1, Nomega = 1000,
                  Nper=1, gamma = 0.0001, Nk = 0, dEk = 0.001, Ek = -1):
        '''Run hdiag8 program in path "path" to diagonalize the basis for a
        specific electric field. In addition, the program calculates:
        Gain with Fermi's golden rule in re-diagonalized Wannier-Stark basis
        (from omega0 to omegaf eV with Nomega points, Nper periods, and a FWHM
        of gamma eV);
        k-resolved occupation function fk for each level in Wannier-Stark
        (with Nk points separated by dEk eV);
        Spectral function A_i(E,E_k) and Im(G^<_i(E,E_k) in Wannier-Stark
        basis (for Ek eV).
        '''

        if omega0 == -1:
            commands = str(zshift)+"\n-1\n0\n-1\n"
        else:
            commands = str(zshift)+"\n"+str(omega0)+"\n"+str(omegaf)+"\n"+str(Nomega)
            commands += "\n"+str(Nper) +"\n"+str(gamma)+"\n"+str(Nk)+"\n"

        if Nk > 0:
            commands += str(dEk) + "\n"
        commands += str(Ek)

        proc = su.dispatch(self.proghdiag, [], path, infile = subprocess.PIPE, outfile = subprocess.PIPE, errfile = subprocess.PIPE)
        out, err = proc.communicate(commands)
        dbg.debug("from hdiag: " + str(out) + str(", ") + str(err), dbg.verb_modes["chatty"], self.__class__)
        su.waitforproc(proc, 0.1)
        return proc

    def checkWSdens(self, einspath, tolerance = 0.1):
        '''Check if doping density in diagonalized basis matches
        the intended doping. Default tolerance is 10%.
        '''

        with open(einspath + "/wslevelsRediag.dat", 'r') as f:
            for line in f:
                dopdiag = -1
                doping = -1
                for i in range( len(line.split()) ):
                    word = line.split()[i]
                    if word == 'indices=':
                        dopdiag = float( line.split()[i+1] )
                        doping = float( line.split()[i+3] )
                        break
                if dopdiag != -1:
                    break

        if dopdiag == -1:
            return False
        elif (abs(dopdiag - doping) > tolerance*doping):
            return False
        else:
            return True


    def getWSdata(self,path):
        '''Get the Wannier-Stark data from wslevels.dat and matrixWS.dat.'''
        levels = []
        with open(path+"/wslevels.dat", 'r') as ff:
            # read eFd
            line = ff.readline().split()
            efd = line[2]
            # read 3 lines of comments
            for _ in range(0,3):
                next(ff)
            for _ in range(0,self.numpar["Nstates"]):
                line = ff.readline().split()
                data = []
                data.append(path)
                data.append(float(efd))
                for element in line:
                    data.append(float(element))
                levels.append(data)
        dipoles = []
        with open(path + "/matrixWS.dat",'r') as ff:
            line = ff.readline().split()
            line = ff.readline().split()
            nlevels = int(line[0])
            for _ in range(3+nlevels*2):
                next(ff)
            for _ in range(nlevels):
                row = []
                line = ff.readline().split()
                for i in range(nlevels):
                    row.append( float(line[i] ) )
                dipoles.append(row)

        return levels, dipoles

    def runBandplot(self, path, structure, valence = 0, localdata = 1, eFd = 0, plotWS = 1,
                    whichWS = 1, square = 1, renorm = 0.3, pmin = -1, pmax = 1, zmin = None, zmax = None):
        '''Runs the bandplot program and saves a plottable file "bandplot.dat".'''

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
            raise Exception("mpameter not valid: plot local data = "+localdata)
        commands += str(zmin) + "\n" + str(zmax) + "\n" + str(pmin) + "\n" + str(pmax) + "\n" + str(square) + "\n" + str(renorm) + "\n"

        proc = su.dispatch(self.progbandplot, [], path, infile = subprocess.PIPE)
        proc.communicate(commands)
        return proc

    def getMerit(self,structure,path, only_conv = True):
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

        path = os.path.join( path , structure.dirname , self.datpath )

        results = []

        try:
            with open(path+"/negft.dat",'r') as f:
                for line in f:
                    try:
                        if '#' in line or only_conv and (line.split()[Inegf.idat.get("ierror")].endswith('1')\
                                or line.split()[Inegf.idat.get("konv")] == 'NaN'):
                            continue
                    except(IndexError): # missing # in negft.dat, with intel compiler
                        continue
                    results.append(line.split())
                if results == []:
                    return "NO CONV"
            f.closed

        except (OSError, IOError):
            print("\nWARNING: Could not find directory: " + path)
            return "ERROR"
        except ValueError:
            print("\nWarning: Error when getting results from: " + path)
            return "ERROR"

        values = []
        for i in range(0, len(Inegf.idat)):
            datval = []
            for row in results:
                datval.append(float(row[i]))
            values.append(datval)

        if self.merit == Interface.merits.get("max gain") :
            # interpolate results:
            iom = Inegf.idat["omega"]
            igain=Inegf.idat["gain"]
            try:
                if len(values[1]) > 2:
                    if len(values[1]) > 3:
                        kind = "cubic"
                    elif len(values[1]) == 3:
                        kind = "quadratic"
                    x = np.linspace(values[iom][0],values[iom][-1])
                    f = interp1d(values[iom], values[igain], kind=kind)
                    out = np.max(f(x))
                else:
                    out = np.max(values[igain])
            except ValueError:
                return "ERROR"

        elif self.merit == Interface.merits.get("(max gain)/(current density)"):

            gain = values[Inegf.idat.get("gain")]
            j = values[Inegf.idat.get("j")]
            gainj = []
            [gainj.append(gain[i]/j[i]) for i in range(0,len(gain))]
            out = max(gainj)

        elif self.merit == self.merits.get("estimated gain") :
            try:
                dirlist = su.listdirs(os.path.join(path,"eins"))
            except (OSError, IOError):
                print("WARNING: could not find directory: " + os.path.join(path,"eins"))
                return "ERROR"
            maxgain = []
            for folder in dirlist:
                einspath = os.path.join(path,"eins",folder)
                if( self.checkWSdens(einspath) ):
                    with open( os.path.join(einspath,"/gainFGR.dat") ) as f:
                        for line in f:
                            linestr = line.split()
                            try:
                                maxgain.append( float( linestr[1] ) )
                            except( IndexError ):
                                pass

            out = max(maxgain)

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
            currlist = np.abs( values[Inegf.idat.get("j")] )
            imax = np.argmax(currlist)
            jmax = currlist[imax]
            omegamax = values[Inegf.idat.get("omega")][imax]

            return jmax*np.exp(-np.abs(omegamax - omegat)/a)

        elif self.merit == self.merits["gain integral"]:
            # Optimize the area under the gain curve

            gain = values[Inegf.idat.get("gain")]
            omega = values[Inegf.idat.get("omega")]
            efd = values[Inegf.idat.get("eFd")]

            sum = 0.
            sums = []
            efd0 = efd[0]
            domega = omega[1] - omega[0]
            for i in range(len(gain)):
                if efd[i] != efd0:
                    sums.append(sum)
                    sum  = 0.
                    efd0 = efd[i]
                sum = sum + gain[i]*domega
            sums.append(sum)
            # Scaling output to prevent sharp narrow peaks
            return np.max(sums)/np.max(np.abs(gain))

        elif self.merit == self.merits["gain FWHM"]:
            # Opitmize the width of the gain peak

            if self.target is not None:
                losses = self.target[0]
            else:
                losses = 0.

            gain = values[Inegf.idat.get("gain")]
            omega = values[Inegf.idat.get("omega")]
            efd = values[Inegf.idat.get("eFd")]
            efd0 = efd[0]
            FWHMs = []
            omega_efd = []
            gain_efd = []
            for i in range(len(gain)):
                omega_efd.append(omega[i])
                gain_efd.append(gain[i])
                if efd[i] != efd0 or i == len(gain) -1:
                    efd0 = efd[i]
                    # end of eFd period
                    om0 = None
                    om1 = None
                    for j in range(len(omega_efd)-1):
                        if (gain[j]-losses)*(gain[j+1]-losses) < 0:
                            f = interp1d([gain[j], gain[j+1]],
                                         [omega[j], omega[j+1]])
                            om0 = f(losses)
                        if (gain[j]-losses)*(gain[j+1]-losses) < 0 and om0 is not None:
                            f = interp1d([gain[j], gain[j+1]],
                                         [omega[j], omega[j+1]])
                            om1 = f(losses)
                    if om0 is None or om1 is None or np.max(gain_efd) < 0:
                        FWHMs.append(0.)
                    else:
                        FWHMs.append(abs(om1-om0))
                    omega_efd = []
                    gain_efd = []
                    om1, om0 = None, None

                    out = np.max(FWHMs)
        else:
            print("No such merit function!")

        return out

    def writeMaterial(self,material,nametag,dirpath = None):
        '''Writes the material.inp input file.'''

        if dirpath is None:
            path="material.inp"
        else:
            path = os.path.join( dirpath , "material.inp" )
        with open(path, 'w') as f:
            f.write(nametag + "\n")
            f.write(str(material.params["meff"])+" # meff\n")
            f.write(str(material.params["eps0"])+" # eps0\n")
            f.write(str(material.params["epsinf"])+" # epsinf\n")
            f.write(str(material.params["ELO"])+" # ELO\n")
            #f.write(str(np.abs(material.params["ac"]))+" # ac\n")
            f.write(str(material.params["ac"])+" # deformation potential\n")
            f.write(str(material.params["vlong"])+" # vlong\n")
            f.write(str(material.params["massdens"])+" # mass density\n")
            f.write(str(material.params["molV"])+ " # mol volume\n")


    def writeWannier(self,struct,dirpath=None):
        '''Writes the wannier8.inp input file.'''

        if dirpath is None:
            path="wannier8.inp"
        else:
            path = os.path.join(dirpath,"wannier8.inp")

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

    def writeNegftInp(self,pathwann,patheins,dirpath=None, numpar=None):
        '''Writes the negft8.inp input file.'''

        if numpar is None:
            numpar = self.numpar

        if numpar["boolEins"]:
            readeins = ".TRUE."
        else:
            readeins = ".FALSE."
        if dirpath is None:
            path = "negft7.inp"
        else:
            path = dirpath+"/negft7.inp"
        try:
            with open(path,"w") as f:
                f.write(str(numpar["Nper"])+" "+str(numpar["Nstates"])+" "+str(numpar["NE"]) + " " + str(numpar["Nk"]))
                f.write(" # Nper, Nnu, NE, Nk\n")
                f.write(str(numpar["Nh"])+" "+str(numpar["Igauge"]))
                f.write(" #\n")
                f.write(str(numpar["Tlattice"]))
                f.write(" #\n")
                f.write(str(numpar["efield0"])+" "+str(numpar["defield"])+" "+str(numpar["Nefield"])+" ")
                f.write(" #\n")
                f.write(str(numpar["efac0"])+" "+str(numpar["defac"])+" "+str(numpar["Nefac"])+" ")
                f.write(" #\n")
                f.write(str(numpar["omega0"])+" "+str(numpar["domega"])+" "+str(numpar["Nomega"])+" ")
                f.write(" #\n")
                f.write(str(numpar["Emaxadd"])+ " " + str(numpar["Eminsub"])+" #\n")
                f.write(str(numpar["gen"])+ " #\n")
                f.write(str(numpar["maxits"])+" "+str(numpar["Iconv"])+" "+str(numpar["Bei"])+" " \
                        + str(numpar["Nhist"]) + " #\n")
                f.write(pathwann + "/\n")
                f.write(readeins+"\n")
                f.write(patheins+"\n")
                if(numpar["use-poisson"]):
                    f.write(".TRUE.")
                else:
                    f.write(".FALSE.")
                f.write(" #\n")
                if(numpar["boolPrinc"]):
                    f.write(".TRUE.")
                else:
                    f.write(".FALSE.")
                f.write(" #\n")
                if(numpar["use-e-e"]):
                    f.write(".TRUE.")
                else:
                    f.write(".FALSE.")
                f.write(" #\n")

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
            try:
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
            except FileNotFoundError:
                print(f"Warning in Inegf.loadStructures: structure not found in {folder}")
        return structures


    def getresults(self, structure, path, datpath = None):
        """Returns results from the NEGF program.

        Paramters:

        structure: Structure
            The structure for which to return results.
        path: String
            The base path where the structure directory resides.
        datpath: String
            (Optional) The execution path of the NEGF program. Defaults to self.datpath.

        Returns:

        negft: numpy.Array
                Matrix with NEGF results.
        """

        if datpath is None:
            datpath = self.datpath

        if path[-4:-1] == '.dat':
            pass
        elif structure is None:
            path = path + "/negft.dat"
        else:
            path = path + "/" + structure.dirname + "/" + datpath + "/negft.dat"

        negft = []
        with open(path, "r") as file:
            for l in file:
                row = []
                try:
                    row = [ float( l.split()[i] ) for i in range(5) ]
                except (ValueError, IndexError):
                    # These are commented lines
                    continue
                negft.append(row)
        negft = np.array(negft)
        negft = negft[negft[:,0].argsort()]

        return negft



    def plotresults(self, structure, path, label = None):

        negft = self.getresults(structure,path)

        pl.figure('IV')
        pl.plot(negft[:,0]*1000, negft[:,3], '-', label=label)
        pl.xlabel('Bias (mV/period)')
        pl.ylabel('Current density (A/cm^2)')
        pl.legend()
        pl.figure('Gain')
        pl.plot(negft[:,1]*1000, negft[:,4], '-+', label=label)
        pl.xlabel('Energy (meV)')
        pl.ylabel('Gain (1/cm)')
        if label is not None:
            pl.legend()

        return negft


    def plotbands(self, structure, path, einspath = None):
        '''
        Plot the band structure and Wannier-Stark states for a structure.
        Must first have run gatherResults() for structure.

        Parameters

        structure: Structure for which to plot the bands.
        path: The path to the base of the tree
        einspath (optional): specifies which parameters to plot for
        '''

        path = path + "/" + structure.dirname + "/" + self.datpath + "/eins/"

        if einspath == None:
            dirlist = su.listdirs(path)
        else:
            dirlist = [einspath]

        Nplots = len(dirlist)
        Nrows = int( np.sqrt(Nplots) )
        Ncols = int( Nplots/Nrows + 0.5 )

        z_all = []
        Ec_all = []
        wavef_all = []

        iplot = 1

        for dir in dirlist:
            with open(path + "/" + dir + "/bandplot.dat", 'r') as bandplot:

                eFd = float( bandplot.readline().split()[2] )
                for _ in range(2):
                    bandplot.readline()   # read comments
                nnu = int( bandplot.readline().split()[1] )
                nz  = int( bandplot.readline().split()[2] )


                z = []
                Ec = []
                wavef = []
                [wavef.append([]) for _ in range(nnu)]
                for l in bandplot:
                    data = l.split()
                    z.append( float( data[0] ) )
                    Ec.append( float( data[1] ) )
                    for nu in range(2, len(data) ):
                        wavef[nu-2].append( float( data[nu] ) )

                z_all.append(z)
                Ec_all.append(Ec)
                wavef_all.append(wavef)

                if Nplots > 1:
                    pl.subplot(Nrows, Ncols, iplot)
                    iplot += 1
                pl.plot(z, Ec)
                for nu in range(0, nnu):
                    pl.plot(z, wavef[nu] )



        return z_all, Ec_all, wavef_all

    def plotGainFGR(self, structure, path, einspath = None, multiplot = True, plot1D = True, plot2D = True,
                   cmap='jet',vmin=None, vmax=None, yrange=None, xrange=None, ydata = 'eFd'):
        '''
        Plot the gain calculated from Fermi's golden rule for a structure.
        Must first have run gatherResults() for structure.

        Parameters

        structure: Structure
            Structure for which to plot the bands.
        path: string
            The path to the base of the tree
        einspath: string
            specifies which parameters to plot for
        multiplot:Boolean
            Make one subplot per eins path
        plot1D: Boolean
            Make a 1D plot with overlain gain curves
        plot2D: Boolean
            Make a 2D color map
        cmap: string
            Which cmap to use in 2D plot
        vmin: float
            Minimun limit on color axis
        vmax: float
            Maximum limit on color axis
        yrange: Tuple
            Range of values on y-axis (eFd or current)
        xrange: Tuple
            Range of values on x-axis (frequency)
        ydata: string
            ydata = 'eFd' gives bias on y-axis, and ydata = 'Current' gives current density.

        Returns:

        om_all: array[array[float]]
            All frequency value arrays for all eins directories.

        g_all: array[array[float]]
            All gain value arrays for all eins directories.
        '''

        path = path + "/" + structure.dirname + "/" + self.datpath + "/eins/"

        if einspath == None:
            dirlist = su.listdirs(path)
        else:
            dirlist = [einspath]

        Nplots = len(dirlist)
        Nrows = int( np.sqrt(Nplots) )
        Ncols = int( Nplots/Nrows + 0.5 )

        om_all = []
        g_all = []

        iplot = 1

        efd_all = []


        if ydata == 'Current':
            negft = []
            with open(path + "/../negft.dat", "r") as file:
                for l in file:
                    row = []
                    try:
                        row = [ float( l.split()[i] ) for i in range(5) ]
                    except (ValueError, IndexError):
                        # These are commented lines
                        continue
                    negft.append(row)
            negft = np.array(negft)
            negft = negft[negft[:,0].argsort()]

        for dir in dirlist:
            try:
                with open(path + "/" + dir + "/gainFGR.dat", 'r') as gain:
                    efd_all.append(float(dir.split(".")[0]))

                    gain.readline()   # one comment line

                    om = []
                    g = []
                    for l in gain:
                        data = l.split()
                        om.append( float( data[0] ) )
                        g.append( float( data[1] ) )

                    om_all.append(om)
                    g_all.append(g)
            except FileNotFoundError:
                pass

        # sort data according to bias
        index_sort = np.argsort(efd_all)
        efd_all = np.array(efd_all)[index_sort]
        g_all = np.array(g_all)[index_sort]
        om_all = np.array(om_all)[index_sort]

        if ydata == 'Current':
            j_all = []
            for i in range(len(efd_all)):
                for j in range(len(negft)):
                    if efd_all[i] == negft[j][0]*1000:
                        j_all.append(negft[j][3])
                        break


        if plot1D:

            pl.figure('FGR_1D')
            for i in range(len(efd_all)):
                # make 1D plots
                pl.figure('FGR_1D')
                if Nplots > 1 and multiplot:
                    pl.subplot(Nrows, Ncols, iplot)
                    iplot += 1
                pl.plot(om_all[i], g_all[i], label = str(efd[i]))
                pl.xlabel('Energy (meV)')
                pl.ylabel('Gain (1/cm)')

        if plot2D:
            pl.figure('FGR_2D')
            # make 2D plot bias vs. frequency (assumes all freq. arrays are same length
            pl.figure('FGR_2D')
            if ydata == 'Current':
                ploty = j_all
            elif ydata == 'eFd':
                ploty = efd_all
            pl.pcolormesh(om_all[0]*1000, ploty, g_all, vmin = vmin, vmax = vmax, cmap = cmap)
            pl.xlabel('Energy (meV)')
            if ydata == 'eFd':
                pl.ylabel('Bias (mV/period)')
            elif ydata == 'Current':
                pl.ylabel('Current density (A/cm^2)')
            if(xrange is not None):
                pl.xrange(xrange)
            if(yrange is not None):
                pl.yrange(yrange)
            cbar = pl.colorbar()
            cbar.set_label('Gain (1/cm)')


        return om_all, g_all, efd_all

    def plotResolve(self, path, structure = None, einspath = None,
                    plotak = True, plotbands = False, vmax_dens = None,
                    vmin_dens = None, vmin_curr = None, vmax_curr = None,
                    cmap_dens = 'hot', cmap_ak = 'cool', WScolor = 'b',
                    colorbar = True):

        if structure is not None:
            path = path + "/" + structure.dirname

        path = path + "/" + self.datpath + "/eins/"

        if einspath == None:
            dirlist = su.listdirs(path)
            subdirs = []
            for d in dirlist:
                if "resolve8.info" in su.ls(path + "/" + d):
                    subdirs.append(d)
            dirlist = subdirs
        else:
            dirlist = [einspath]

        Nplots = len(dirlist)
        Nrows = int( np.sqrt(Nplots) )
        Ncols = int( Nplots/Nrows + 0.5 )

        iplot = 1

        fig1 = pl.figure() # Figure for current
        fig2 = pl.figure() # Figure for density

        for dir in dirlist:
            pathdir = path + "/" + dir
            bandplot = np.transpose( np.loadtxt(pathdir + "/bandplot.dat") )
            if plotak:
                try:
                    ak0 =  np.loadtxt( pathdir + "/resolve_ak.dat")
                except FileNotFoundError:
                    plotak = False
            dens = np.loadtxt( pathdir + "/resolve_dens.dat")
            curr = np.loadtxt( pathdir + "/resolve_curr.dat")

            info = open(pathdir + "/resolve8.info","r")

            for line in info:
                if "columns" in line.split():
                    Nz, zmin, zmax = int(line.split()[0]), float(line.split()[6]), float(line.split()[8])
                elif "rows" in line.split():
                    NE, Emin, Emax = int(line.split()[0]), float(line.split()[6]), float(line.split()[8])

            z = np.linspace(zmin,zmax,Nz)
            E = np.linspace(Emin*1000,Emax*1000,NE)

            N = 512 # number of colors to use
            zeros = [[0,0],[0,0]]

            pl.figure(fig1.number)
            if Nplots > 1:
                pl.subplot(Nrows, Ncols, iplot)
                #iplot += 1
            pl.title(dir)
            pl.plot(bandplot[0],bandplot[1],WScolor)
            if plotak:
                pl.contour(z,E,ak0,cmap=cmap_ak)
            if plotbands:
                for i in range(2,len(bandplot)):
                    pl.plot(bandplot[0],bandplot[i],WScolor)

            #f=pl.contourf(z,E,curr,N, cmap = 'hot', vmin = vmin, vmax = vmax)
            pl.pcolormesh(z, E, curr, cmap = cmap_dens, vmin = vmin_curr, vmax = vmax_curr)
            if colorbar: pl.colorbar()
            pl.xlim(zmin,zmax)
            pl.xlabel("z (nm)")
            pl.ylabel("E (meV)")
            Emin2 = pl.ylim()[0]
            #pl.contourf([zmin,zmax],[Emin2,Emin*1000],zeros,cmap='hot',levels=f.levels)
            pl.pcolormesh([zmin,zmax],[Emin2,Emin*1000],zeros,cmap=cmap_dens, vmin = vmin_curr, vmax = vmax_curr)

            pl.figure(fig2.number)
            if Nplots > 1:
                pl.subplot(Nrows, Ncols, iplot)
                iplot += 1
            pl.title(dir)
            pl.plot(bandplot[0],bandplot[1],WScolor)
            if plotak:
                pl.contour(z,E,ak0,cmap=cmap_ak)
            if plotbands:
                for i in range(2,len(bandplot)):
                    pl.plot(bandplot[0],bandplot[i],WScolor)
            #f=pl.contourf(z,E,dens,N, cmap = 'hot')
            pl.pcolormesh(z,E,dens, cmap = cmap_dens, vmin=vmin_dens, vmax=vmax_dens)
            if colorbar: pl.colorbar()
            pl.xlim(zmin,zmax)
            pl.xlabel("z (nm)")
            pl.ylabel("E (meV)")
            Emin2 = pl.ylim()[0]
            #pl.contourf([zmin,zmax],[Emin2,Emin*1000],zeros,cmap='hot',levels=f.levels)
            pl.pcolormesh([zmin,zmax],[Emin2,Emin*1000],zeros,cmap=cmap_dens,vmin=vmin_dens, vmax=vmax_dens)

        return fig1, fig2


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
        elevel = []
        nlevel = []
        glevel = []
        # TODO: only takes the first eins folder
        for level in structure.wslevels[0]:
            elevel.append( level[4] )
            nlevel.append( level[5] )
            if gamma is None:
                glevel.append( level[6] )
            else:
                glevel.append(gamma)

        dipoles = structure.dipoles[0]
        dopdens = 0.
        for d in structure.dopings:
            dopdens += (d[1]-d[0])*d[2]
        dopdens /= structure.length

        if len(elevel) < 3 or len(dipoles) < 3 or self.checkWSdens(structure.wslevels[0][0][0], 0.1) == False:
            return -1

        if E1 is None:
            E1 = elevel[2] - elevel[0]
        if E2 is None:
            E2 = elevel[1] - elevel[0]

        E3 = E1-E2


        chi2 =  (nlevel[0]-nlevel[2])/( E1-elevel[2]+elevel[0] -
                    (glevel[0] + glevel[2])*1j/2.)
        chi2 += (nlevel[0]-nlevel[1])/(-E2+elevel[1]-elevel[0] -
                    (glevel[0] + glevel[1])*1j/2.)
        chi2 *= dipoles[0][1]*dipoles[0][2]*dipoles[2][1]
        chi2 *= structure.length**3 # from Z0/d -> nm
        chi2 /= (E3-elevel[2]+elevel[1] -
                    (glevel[1] + glevel[2])*1j/2.)
        chi2 *= -const.qe/const.eps0/structure.length*1e4*1e-27*1e9*1e12

        self.chi2 = chi2

        return np.abs( chi2 )
