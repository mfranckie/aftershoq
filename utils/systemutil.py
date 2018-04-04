'''
Created on 26 Jan 2018

@author: martin
'''

from subprocess import call, Popen, PIPE
import os
import time
from utils.debug import Debugger as dbg
from tempfile import TemporaryFile as tmp

class SystemUtil(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
        
    @staticmethod
    def mkdir(name):
        call(["mkdir",name])
        
    @staticmethod
    def rmdir(name):
        call("rm -r " + name,shell=True)
        
    @staticmethod
    def ls(name):
        call("ls -l "+name,shell=True)
    
    @staticmethod
    def listdirs(path):
        # find . -type d -depth 1: lists directories only
        output = []
        [output.append(d) for d in os.listdir(path) if os.path.isdir(path+"/"+d)]
        
        return output
    
    @staticmethod
    def waitforproc(proc,delay = 0.1):
        while proc.poll() == None:
            time.sleep()
        
    @staticmethod
    def dispatch(prog,args,dirpath=None, infile = None, outfile = None, errfile = None):
        
        if infile is None:
            infile = tmp()
        if outfile is None:
            outfile = tmp()
        if errfile is None:
            errfile = tmp()
        
        progargs=[prog]
        
        for a in args:
            progargs.append(a)
        
        if dirpath is None:
            dirpath = "./"
        
        dbg.debug( "<<<< Dispatching program: " + str(progargs) + " from " + dirpath,
                   dbg.verb_modes["chatty"], SystemUtil)
        
        process=Popen(progargs,cwd=dirpath,close_fds=True,stdout=outfile,stderr=errfile,stdin=infile)
        
        dbg.debug(  " with pid="+str(process.pid)+" >>>>\n" ,
                    dbg.verb_modes["chatty"])
        dbg.flush()
        return process
