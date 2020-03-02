'''
Created on 18 Jun 2018

@author: Martin Franckie

'''

from aftershoq.numerics import OptimizerND
import numpy as np
from matplotlib import pyplot as pl
import GPyOpt
import aftershoq.utils.debug as dbg

class MDGaussopt(OptimizerND):
    '''
    OptimizerND class for multi-dimensional Gaussian process optimization.
    After initializastion the optimal points to evaluate next are given by
    nextstep(). A fully automatic minimization can be performed by calling
    minimize().
    '''

    def __init__(self, tolerance, maxiter, procmax, limits, x0 = [], y0 = [],
                 **argsGPy):
        '''Constructor.

        Parameters

        tolerance: the absolute tolerance as a distance along the HilbertCurve
        maxiter: the maximum iterations for the optimization.
        procmax: the maximum points to be evaluated concurrently.
        limits: list of the domain limits for each parameter.
        x0 (optional): Initial training data parameters
        y0 (optional): Initial training data function values
        argsGPy (optional): Arguments to pass to GPyOpt. If None, standard
        arguments for parallel optimization are used.
        '''

        self.xmin = 0
        self.iter = 0

        if len(argsGPy)==0:
            argsGPy = {'acquisition_type':'EI','normalize_Y':False,'batch_size':procmax,'num_cores':procmax,
              'maximize':False,'evaluator_type':'local_penalization'}

        domain = []
        for i in range(len(limits)):
            domain.append({'name': 'var'+str(i), 'type': 'continuous', 'domain': limits[i]})

        argsGPy.update({'domain':domain})

        print(f'Creating GP optimizer with parameters:')
        print(f'{argsGPy}')

        self.argsGPy = argsGPy

        # Put x and y values on correct form
        x0 = np.array(x0)
        y0 = np.array(y0)

        if y0.ndim != 2:
            y0 = np.reshape(y0,(len(y0),1))

        self.GP_model = GPyOpt.methods.BayesianOptimization(f = None,
                            X = x0, Y = y0, **argsGPy)

        if len( x0 ) > 0:
            self.xmax = np.max(x0)
        else:
            self.xmax = 0

        OptimizerND.__init__(self, tolerance, maxiter, procmax, limits, x0, y0)

    def addpoints(self, newx, newy):

        #if ( type(self.x) is np.ndarray):
        #    self.x = np.ndarray.tolist(np.squeeze( self.x ))
        #    self.y = np.ndarray.tolist(np.squeeze( self.y ))

        OptimizerND.addpoints(self, newx, newy)

        if self.y.ndim != 2:
            self.y = np.reshape(self.y,(len(self.y),1))

        self.GP_model = GPyOpt.methods.BayesianOptimization(f = None,
                            X = self.x, Y = self.y, **self.argsGPy)


    def nextstep(self):

        newx = self.GP_model.suggest_next_locations()

        # Check convergence:
        self.check_conv()

        return newx

    def minimize(self, model, sgenerator, pathwd, pathresults = None, seq=None):

        if pathresults is None:
            pathresults = pathwd

        dbg.debug("Starting minimization\n", dbg.verb_modes["verbose"],self)
        niter = 0
        while self.converged == 0:
            niter += 1
            dbg.debug("Iteration " + str(niter) + "\n",
                      dbg.verb_modes["verbose"],self)

            newx = self.nextstep()
            sgenerator.generateStructures(newx)
            if seq is not None:
                t = model.runStructSeq(sgenerator.structures[-len(newx):], pathwd,seq)
                for tt in t: tt.join()
            else:
                model.runStructures(sgenerator.structures[-len(newx):], pathwd)
                model.waitforproc(60)
            newy = []
            newx_success = []
            xi = 0
            model.gatherResults(sgenerator.structures[-len(newx):], pathwd, pathresults = pathresults)
                
            self.addEvaldPoints(model, sgenerator, pathwd, sgenerator.structures[-len(newx):])
                
            #self.writeresults(pathresults, "mdgauss.log")

        dbg.debug("Minimization finished with convergence: " + str(self.converged) + "\n",
                  dbg.verb_modes["verbose"],self)
        dbg.flush()
        return self.converged



    def getbest(self):
        index = np.argmin(self.y)
        return [self.x[index], self.y[index] ]



    def writeresults(self, pathresults, filename):
        out = []
        for i in range(len(self.x)):
            row = np.array( [ self.y[i], self.x[i] ] )
            out.append(row)
        np.savetxt(os.path.join(pathresults,filename),out)
