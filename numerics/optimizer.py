'''
Created on 18 Jun 2018

@author: martin
'''

import numpy as np

class Optimizer1D(object):
    '''
    classdocs
    '''


    def __init__(self, tolerance, maxiter, procmax, x0 = [], y0 = []):
        '''
        tolerance: optimization will converge after Delta x < tolerance.
        r: parameter of the optimization scheme.
        maxiter: optimization will terminate after maxiter iterations.
        procmax: maximum number of processors that can be used. This
                is the number of x-values that will be provided by the
                algorithm
        x0: Initial x-values that are already evaluated
        y0: Initial results of the evaluations at x0
        '''
        self.pmax = procmax
        self.x = []
        self.y = []
        self.tol = tolerance
        self.maxits = maxiter
        self.t = -1
        self.iter = 0
        self.converged = 0
        if(len(x0) > 0 and len(y0) > 0):
            self.addpoints(x0, y0)
        
    def addpoints(self,newx,newy):
        '''
        Add the newly evaluated y-values newy at the coordinates newx.
        The new results will be sorted and convergence checked.
        '''
        
        [self.x.append(xx) for xx in newx]
        [self.y.append(yy) for yy in newy]
        
        self.t = np.argmin(self.y)
        self.check_conv()
    
    def nextstep(self):
        pass
    
    def check_conv(self):
        '''
        Check if the minimization has converged.
        Returns:
        0 if it has not converged yet
        1 if it has converged
        -1 if it has reached the maximum allowed iterations
        '''
        
        if self.iter > self.maxits:
            self.converged = -1
            return
        if self.t == len(self.x)-1:
            start = self.t-2
        else:
            start = self.t-1
        for i in range(start,start+2):
            conv = abs(self.x[i]-self.x[i+1])
            if conv < self.tol:
                self.converged = 1
                return
        self.converged = 0
        
    def minimize(self, model, sgenerator, pathwd, pathresults = None):
        pass
    
    def minimize_parameters(self, model, hutil):
        pass
    
    def getbest(self):
        pass