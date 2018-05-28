'''
Created on 14 May 2018

@author: martin
'''
from interface import Interface
import numpy as np
import random as rm

class NDtestfunc(Interface):
    
    '''
    ND = number of dimensions for parameter space
    '''
    def __init__(self, ND, C = None):
        
        self.ND = ND
        self.NC = 4
        if C is None:
            C = []
            for _ in range (0, self.NC):
                A = []
                for i in range (0, ND+1):
                    A.append([])
                    [A[i].append(rm.random()*2-1) for _ in range(0,ND+1)]
                C.append(A)
        self.C = C
        
    '''
    Returns the objective evaluated for an N-dimensional parameter list
    @x: N-dimensional list of parameters.
    '''
    def testfunc(self, x):
        fac = 0.30
        a = b = 0
        for i in range(0,self.ND+1):
            for j in range(0,self.ND+1):
                ai = 1; bi = 1
                for n in range(0, self.ND):
                    ai *= np.sin((n+i)*np.pi*x[n]*fac)
                    bi *= np.cos((n+j)*np.pi*x[n]*fac)
                A = self.C[0]; B = self.C[1]
                C = self.C[2]; D = self.C[3]
                a += A[i][j]*ai + B[i][j]*bi
                b += C[i][j]*ai + D[i][j]*bi
        return -np.sqrt(a**2 + b**2)
    
    def runStructures(self, structures, numpar, path):
        pass
    
    def getMerit(self, params):
        return self.testfunc( params )
    