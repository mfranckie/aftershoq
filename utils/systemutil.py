'''
Created on 26 Jan 2018

@author: Martin Franckie

Module for handling interactions with the system, such as running external
programs.
'''

from subprocess import call, Popen, PIPE
import os
import time
import utils.debug as dbg
from tempfile import TemporaryFile as tmp

def mkdir(name):
    '''Make directory "name", non-recursively. Calls mkdir.'''
    
    call(["mkdir",name])
    
def rmdir(name):
    '''Delete directory "name". Calls rm.'''
    call("rm -r " + name,shell=True)
    
def ls(path):
    '''List files and directories in path. Calls ls.'''
    call("ls -l "+name,shell=True)

def listdirs(path):
    '''List directories only in path.'''
    output = []
    [output.append(d) for d in os.listdir(path) if os.path.isdir(path+"/"+d)]
    
    return output

def waitforproc(proc,delay = 0.1):
    '''Wait for the process proc to finish. Sleep with delay seconds.'''
    while proc.poll() == None:
        time.sleep()
    
def dispatch(prog,args,dirpath=None, infile = None, outfile = None, errfile = None):
    '''Dispatch program prog with arguments args.
    
    Optional parameters:
    dirpath: Path to working directory.
    infile:  File with inputs to the program. Use subprocess.PIPE for input stream.
    outfile: File for program output. Use subprocess.PIPE for output stream.
    errfile: File for program output. Use subprocess.PIPE for output stream.
    
    Returns the Popen() process.
    
    Note that if all outputs are piped, the system can cause the program to crash
    since the pipes might not be properly closed between calls.
    '''
    
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
               dbg.verb_modes["chatty"])
    
    process=Popen(progargs,cwd=dirpath,close_fds=True,stdout=outfile,stderr=errfile,stdin=infile)
    
    dbg.debug(  " with pid="+str(process.pid)+" >>>>\n" ,
                dbg.verb_modes["chatty"])
    dbg.flush()
    return process
