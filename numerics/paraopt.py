'''
Created on 13 Mar 2018

@author: martin
'''

import numpy as np
from optimizer import Optimizer1D
from utils.debug import Debugger as dbg

class Paraopt(Optimizer1D):
    '''
    Class for 1D function minimization with multiple parallel processors.
    After Strongin and Sergeyev "Global Optimization with Non-Convex
    Constraints - Sequencial and Parallel Algorithms", Kluwer Academic
    Publishers, London (2000).
    '''


    def __init__(self, tolerance, r, maxiter, procmax, x0 = [], y0 = []):
        super(Paraopt,self).__init__(tolerance, maxiter, procmax, x0, y0)
        self.r = r
        
    def addpoints(self, newx, newy):
        Optimizer1D.addpoints(self, newx, newy)
        
        # sort results:
        si = np.argsort(self.x)
        self.x.sort()
        y2 = []
        [y2.append(self.y[yi]) for yi in si]
        self.y = y2
        
    def addEvaldPoints(self, model, sg, path, coords):
        '''
        Convert the evaluated points in N-dim. paramtere space in the 
        structure generator "sg", already evaluated with the Inteface 
        "model", to points along the Hilbert curve and add them.
        '''
        # collect results from trial points
        x0 = []
        [x0.append( sg.hutil.interp_dist_from_coords( c ) ) for c in coords]
        x0 = np.array(x0)
        x0.sort()
        x0 = x0.tolist()
        
        y0 = []
        xi = 0
        for i in range(0,len(x0)):
            try:
                y0.append( -float(model.getMerit(sg.structures[i],path)) )
            except( ValueError ):
                del x0[xi]
                xi-=1
            xi +=1
            
        self.addpoints(x0,y0)
        
        return x0,y0
        
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
        if(mu == 0):
            mu = 1
        
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
        
        p = min(k+1,self.pmax)
        
        for i in range(0,p):
            t = si[k-i]+1
            val = (self.x[t]+self.x[t-1])/2
            
            if t > 1 and t <k+1:
                val -= (self.y[t]-self.y[t-1])/2/self.r/mu
            
            xnext.append( val )
        
        self.t = np.argmin(self.y)
        return xnext
    
    def getbest(self):
        '''
        Returns the so-far minimal coordinates.
        '''
        self.t = np.argmin(self.y)
        return [self.x[self.t],self.y[self.t]]
        
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
        
        if pathresults is None:
            pathresults = pathwd
            
        dbg.debug("Starting minimization\n", dbg.verb_modes["verbose"],self)
        niter = 0
        while self.converged == 0:
            niter += 1
            dbg.debug("Iteration " + str(niter) + "\n", 
                      dbg.verb_modes["verbose"],self)
            
            newx = self.nextstep()
            sgenerator.gen_struct_from_hilbert_curve(newx)
            model.runStructures(sgenerator.structures[-len(newx):], pathwd)
            model.waitforproc(0.1)
            newy = []
            xi = 0
            model.gatherResults(sgenerator.structures[-len(newx):], pathwd, pathresults)
            for ss in sgenerator.structures[-len(newx):]:
                try:
                    val = -float(model.getMerit(ss,pathwd))
                    newy.append( val )
                except( ValueError ):
                    del newx[xi]
                    xi -= 1
                xi += 1
        
            self.addpoints(newx,newy)
            self.writeresults(pathresults, "hilbert.log")
        
        dbg.debug("Minimization finished with convergence: " + str(self.converged) + "\n", 
                  dbg.verb_modes["verbose"],self)
        dbg.flush()
        return self.converged
    
    def minimize_parameters(self, model, hutil):
        
        niter = 0
        while self.converged == 0:
            niter += 1
            
            newx = self.nextstep()
            newy = []
            for xx in newx:
                newy.append( -float(model.getMerit(hutil.interp_coords_from_dist(xx)/float(2**hutil.p))) )
        
            self.addpoints(newx,newy)
        
        return self.converged
        
    def writeresults(self, pathresults, filename):
        with open(pathresults + "/" + filename, 'w') as f:
            [f.write( str( self.x[i] ) + " " + str( self.y[i] ) +"\n") for i in range(0,len(self.x))]
        