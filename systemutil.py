'''
Created on 26 Jan 2018

@author: martin
'''

from subprocess import call, Popen, PIPE
import os
import time
from debug import Debugger as dbg

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
    def dispatch(prog,args,dirpath=None):
        
        progargs=[prog]
        
        for a in args:
            progargs.append(a)
        
        if dirpath is None:
            dirpath = "./"
        
        dbg.debug( "<<<< Dispatching program: " + str(progargs) + " from " + dirpath )
        
        process=Popen(progargs,cwd=dirpath,close_fds=True,stdout=PIPE,stderr=PIPE,stdin=PIPE)
        
        dbg.debug(  " with pid="+str(process.pid)+" >>>>\n" )
        dbg.flush()
        return process

    @staticmethod
    def runprogram(prog,args=None,dirpath=None):
        progargs = [prog]
        if dirpath is None:
            dirpath = "./"
        if args is not None:
            [progargs.append(a) for a in args]
            
        dbg.debug( "<<<< Running program: " + str(progargs)  + " from " + dirpath )
        
        p = Popen(progargs,cwd=dirpath,stdout=PIPE,stderr=PIPE,stdin=PIPE)
        
        dbg.debug( " with pid="+str(p.pid)+" >>>>\n" )
        dbg.flush()
        return p
        
