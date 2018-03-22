'''
Created on 13 Mar 2018

@author: martin
'''

import numpy as np
from utils.debug import Debugger as dbg

class Paraopt(object):
    '''
    Class for 1D function minimization with multiple parallel processors.
    After Strongin and Sergeyev "Global Optimization with Non-Convex
    Constraints - Sequencial and Parallel Algorithms", Kluwer Academic
    Publishers, London (2000).
    '''


    def __init__(self, tolerance, r, maxiter, procmax, x0, y0):
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
        self.r = r
        self.maxits = maxiter
        self.t = -1
        self.iter = 0
        self.converged = 0
        self.addpoints(x0, y0)
        
    def nextstep(self):
        '''
        Returns the procmax number of optimal x-values to be evalueated next.
        '''
        self.iter+=1

        xnext = []
        mulist = []
        R = []
        k = len(self.x)-2
        for j in range(1,k+1):
            if (self.x[j]-self.x[j-1]==0):
                mulist = [1]
                break
            else:
                mulist.append(abs(self.y[j]-self.y[j-1])/(self.x[j]-self.x[j-1]))
        mu = max(mulist)
        
        R.append( 2*(self.x[1]-self.x[0]) - 4*self.y[1]/self.r/mu )
        for j in range(2,k+1):
            if (self.x[j]-self.x[j-1]==0):
                R.append( -2*(self.y[j]+self.y[j-1])/self.r/mu )
            else:
                R.append( self.x[j]-self.x[j-1]+abs(self.y[j]-self.y[j-1])**2/mu**2/(self.x[j]-self.x[j-1]) - 2*(self.y[j]+self.y[j-1])/self.r/mu )
                
        R.append( 2*(self.x[k+1]-self.x[k]) - 4*self.y[k]/self.r/mu )
        
        # find t
        # R starts from 0, not 1!
        si = np.argsort(R)
        t = np.argmax(R)+1
        
        p = min(k+1,self.pmax)
        
        for i in range(0,p):
            t = si[k-i]+1
            val = (self.x[t]+self.x[t-1])/2
            
            if t > 1 and t <k+1:
                val -= (self.y[t]-self.y[t-1])/2/self.r/mu
            
            xnext.append( val )
        
        self.t = t
        return xnext
    
    def getbest(self):
        '''
        Returns the so-far minimal coordinates.
        '''
        return [self.x[self.t],self.y[self.t]]
    
    def addpoints(self,newx,newy):
        '''
        Add the newly evaluated y-values newy at the coordinates newx.
        The new results will be sorted and convergence checked.
        '''
        
        [self.x.append(xx) for xx in newx]
        [self.y.append(yy) for yy in newy]
        
        # sort results:
        si = np.argsort(self.x)
        self.x.sort()
        y2 = []
        [y2.append(self.y[yi]) for yi in si]
        self.y = y2
        self.t = np.argmin(self.y)
        self.check_conv()
        
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
    
    
    def minimize(self, model, sgenerator, path):
        
        dbg.debug("Starting minimization\n", dbg.verb_modes["verbose"],self)
        niter = 0
        while self.converged == 0:
            niter += 1
            dbg.debug("Iteration " + str(niter) + "\n", 
                      dbg.verb_modes["verbose"],self)
            
            newx = self.nextstep()
            sgenerator.gen_struct_from_hilbert_curve(newx)
            model.runStructures(sgenerator.structures[-len(newx):], path)
            model.waitforproc(0.1)
            newy = []
            xi = 0
            model.gatherResults(sgenerator.structures,path)
            for ss in sgenerator.structures[-len(newx):]:
                try:
                    val = -float(model.getMerit(ss,path))
                    newy.append( val )
                except( ValueError ):
                    del newx[xi]
                    xi -= 1
                xi += 1
        
            self.addpoints(newx,newy)
        
        dbg.debug("Minimization finished with convergence: " + str(self.converged) + "\n", 
                  dbg.verb_modes["verbose"],self)
        dbg.flush()
        return self.converged
        