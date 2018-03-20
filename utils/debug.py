'''
Created on 16 Mar 2018

@author: martin
'''

class Debugger(object):
    '''
    Debug class, for writing debugging information.
    '''
    verb_modes = {"silent" : 0, "verbose" : 1, "chatty" : 2}
    file = []
    filename='debug.log'

    def __init__(self, verbosity, outfile = None):
        '''
        
        '''
        self.filename = outfile
        self.verbosity = verbosity
    
    @classmethod
    def open(cls, verbosity, outfile = None):
        cls.verbosity = verbosity
        cls.filename = outfile
        if verbosity == cls.verb_modes["silent"] or outfile is None:
            return
        cls.file = open(outfile,'a')
        
    @classmethod
    def close(cls):
        if cls.verbosity == cls.verb_modes["silent"] or cls.file==[]:
            return
        cls.file.close()
    
    @ classmethod
    def debug(cls, message, level = 2, callclass = None):
        if callclass is not None:
            message = str(callclass) + ": " + message
        if cls.verbosity < level:
            return
        elif cls.file == []:
            print message
        else:
            cls.file.write(message)
    
    @classmethod
    def flush(cls):
        cls.close()
        cls.open(cls.verbosity,cls.filename)