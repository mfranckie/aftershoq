'''
Created on 2 Feb 2018

@author: martin
'''

import aftershoq.utils.systemutil as su
import time
import aftershoq.utils.debug as dbg
import subprocess
import os
import psutil

class Platform(object):
    '''
    Run platform base class. All derived classes are to set
    the self.paral member to one of the self.paral_modes:
    "MPI"
    "OMP"
    "SERIAL" (default)
    '''
    
    paral_modes = {"MPI": 0,
             "OMP": 1,
             "SERIAL": 2}


    def __init__(self, name=None):
        '''
        Constructor
        '''
        
        self.paral = self.paral_modes["SERIAL"]
        
    def submitjob(self,prog,args,dirpath,Nproc=None,wtime=None):
        pass
    
    def jobstatus(self,proc):
        pass
    
class Euler(Platform):
    
    
    
    def __init__(self, Nproc, wtime, paral_in = None):
        '''
        Constructor
        '''
        super(Euler,self).__init__("Euler")
        
        if paral_in is None:
            paral_in = self.paral_modes.get("MPI")
        
        self.paral = paral_in
        if self.paral == self.paral_modes.get("OMP"):
            os.environ["OMP_NUM_THREADS"] = str(Nproc)
        self.Nproc = Nproc
        self.wtime = wtime
        self.subcommand = "bsub"
        self.procinfo = []
    
    def submitjob(self, prog, args,dirpath,Nproc=None,wtime=None):
        if Nproc is None:
            Nproc = self.Nproc
        if wtime is None:
            wtime = self.wtime
        
        progargs = []
        jobname = str(dirpath)
        
        if Nproc>1 and self.paral == self.paral_modes.get("OMP"):
            progargs.append("-R")
            progargs.append('"span[ptile='+str(Nproc)+']"')
        
        progargs.append("-n")
        progargs.append(str(Nproc))
        progargs.append("-W")
        progargs.append(str(wtime))
        progargs.append("-J")
        progargs.append(jobname)
        if Nproc>1 and self.paral == paral_modes.get("MPI"):
            progargs.append("mpirun")
        progargs.append(prog)
        [progargs.append(a) for a in args]
        proc = su.dispatch(self.subcommand, progargs, dirpath)
        self.procinfo.append([proc.pid,jobname])
        return proc
    
    def jobstatus(self, proc):
        for info in self.procinfo:
            if info[0]==proc.pid:
                jobname = info[1]
        p = su.dispatch("bjobs",["-J",jobname],errfile=subprocess.PIPE)
        out = p.communicate()
        
        if ("not found" in out[1]) and (proc.poll()!=None):
            dbg.debug( "process ended!" )
            return False
        else:
            return True
        
class Local(Platform):
    """
    Run platform for running local simulations, on a single node or computer.
    
    Parameters:
    Nproc (optional) : If Nproc > 1 is provided, OpenMP will be used.
                        Defaults to 1. 
    """
    
    def __init__(self, Nproc = 1):
        super(Local,self).__init__("Local")
        
        self.Nproc = Nproc
        if Nproc > 1:
            self.paral = self.paral_modes["OMP"]
            os.environ["OMP_NUM_THREADS"] = str(Nproc)
        
        
    def submitjob(self, prog, args, dirpath, Nproc = None, wtime = None):
        proc = su.dispatch(prog, args, dirpath)
        return proc
    
    def submitandwait(self, prog,  args, dirpath):
        proc = su.dispatch(prog, args, dirpath)
        while proc.poll() == None:
            time.sleep(0.1)
        return proc
    
    def jobstatus(self, proc):
        if proc.poll()==None:
            return True
        else:
            return False
        
class MPI(Platform):
    
    def __init__(self, logical = False):
        super(MPI, self).__init__("MPI")
        self.paral = self.paral_modes["MPI"]
        self.Nproc = psutil.cpu_count(logical=logical)
        
    def submitjob(self, prog, args, dirpath, Nproc = None, wtime = None):
        if Nproc is None:
            Nproc = self.Nproc
        
        progargs = []
        progargs.append("-np")
        progargs.append(str(Nproc))
        progargs.append(prog)
        [progargs.append(a) for a in args]
        
        proc = su.dispatch("mpirun",progargs,dirpath)
        
        return proc
        
    def submitandwait(self, prog,  args, dirpath, Nproc = None):
        proc = self.submitjob(prog, args, dirpath, Nproc)
        while proc.poll() == None:
            time.sleep(0.1)
        return proc
    
    
    def jobstatus(self, proc):
        if proc.poll()==None:
            return True
        else:
            return False