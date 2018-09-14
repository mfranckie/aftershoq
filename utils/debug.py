'''
Created on 16 Mar 2018

@author: Martin Franckie

Debug class, for writing debugging information. The level of detail in the
    output is given by the dictionary "verb_modes", containing levels
    "silent", "verbose", and "chatty", in order of increasing detail.
    Output is written to stdout in verbose mode by default, but can be 
    changed by calling "open()".
    
    Write output via the "debug()" method, and call "flush()" to push the
    output to the stream.
'''

import builtins

verb_modes = {"silent" : 0, "verbose" : 1, "chatty" : 2}
file = None
filename='debug.log'
verbosity = verb_modes["verbose"]


def open(verb, outfile = None):
    '''Open a new file with name "outfile" (optional, defualts to ./debug.log)
    with the given verbosity level (from the dictionary "verb_modes").
    '''
    
    global file, filename, verbosity
    
    verbosity = verb
    filename = outfile
    if verbosity == verb_modes["silent"] or outfile is None:
        return
    file = builtins.open(outfile,'a')
    
def close():
    '''Closes the debug file.'''
    
    global file, verbosity
    
    if verbosity == verb_modes["silent"] or file==[]:
        return
    if file is not None:
        file.close()

def debug(message, level = 2, callclass = None):
    '''Prints a debug message.
    message: string with message
    level:   verbosity level from dictionary verb_modes, defualts to "verbose"
    callclass (optional): class from which debut() is called.
    '''
    
    global file, verbosity
    
    if callclass is not None:
        message = str(callclass) + ": " + message
    if verbosity < level:
        return
    elif file is None:
        print(message)
    else:
        file.write(message)

def flush():
    '''Close and re-open debug file to force the output through the stream.'''
    
    global filename, verbosity
    
    close()
    open(verbosity,filename)