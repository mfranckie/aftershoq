'''
Created on 14 May 2018

@author: Martin Franckie
'''
from aftershoq.interface import Interface
import numpy as np
import random as rm
from matplotlib import pyplot as pl

class NDtestfunc(Interface):

    '''Fast functions for testing optimization schemes.'''

    def __init__(self, ND, C = None, function = 'strongin'):
        '''Constructor.

        Parameters

        ND: Number of dimensions of parameter space.
        C (optional): List of parameters defining the functions.
                      Defaults to random numbers.
        '''

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
        self.testfunc = getattr(self, function)

    def strongin(self, x):
        '''
        Returns the objective scalar evaluated for an N-dimensional parameter
        list. Function from Strongin and Sergeev.

        Parameters

        x: N-dimensional list of parameters.
        '''
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

    def rastrigin(self, X, A = 10):
        """Rastrigin ND function.

        Parameters
        ----------
        X : array[float]
            Coordinate to evaluate at.
        A : float
            Function parameter.

        Returns
        -------
        float
            The Rastrigin function value at X.
        """

        return A*len(X) + sum([(x**2 - A * np.cos(2 * np.pi * x)) for x in X])

    def runStructures(self, structures, numpar, path):
        pass

    def getMerit(self, params):
        '''
        Returns the objective scalar evaluated for an N-dimensional parameter
        list.

        Parameters

        x: N-dimensional list of parameters.
        '''
        return self.testfunc( params )

    def plt(self, x = None, y = None, cmin= None, cmax=None, cmap = 'hot'):

        if x is None:
            x = np.linspace(0,1)
        if y is None and self.ND > 1:
            y = np.linspace(0,1)

        z = []

        if self.ND > 1:
            coord = []
            for i in range(self.ND):
                coord.append(0.)
            for xx in x:
                row = []
                coord[1] = xx
                for yy in y:
                    coord[0] = yy
                    row.append(self.testfunc(coord))
                z.append(row)

        #pl.contour(x,y,z)
        pl.pcolormesh(x,y,z,vmin=cmin, vmax=cmax, cmap=cmap)
        cbar = pl.colorbar()
        return cbar
