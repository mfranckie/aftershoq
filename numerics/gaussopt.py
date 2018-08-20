'''
Created on 18 Jun 2018

@author: Martin Franckie

'''

from optimizer import Optimizer1D
import numpy as np
from matplotlib import pyplot as pl
import scipy.optimize as so
import utils.debug as dbg

class Gaussopt(Optimizer1D):
    '''
    Optimizer1D class for Gaussian process optimization. After initializastion
    the optimal points to evaluate next are given by nextstep(). A fully
    automatic minimization can be performed by calling minimize() or 
    minimize_parameters().
    '''

    def __init__(self, tolerance, maxiter, procmax, x0 = [], y0 = [], sigma = 4, l = 2, sigma_noise = 0., padding = 10.):
        '''Constructor.
        
        Parameters
        
        tolerance: the absolute tolerance as a distance along the HilbertCurve
        maxiter: the maximum iterations for the optimization.
        procmax: the maximum points to be evaluated concurrently.
        x0 (optional): Initial training data parameters
        y0 (optional): Initial training data function values
        sigma: standard deviation of GP (default = 4).
        l: correlation length of GP (default = 2).
        sigma_noise: noise level of GP (default = 0). If 0, noise will not be
        updated.
        padding: minimal distance for evaluating two concurrent points along 
        HilbertCurve.
        '''
        
        self.K = [] # covariance matrix
        self.Kt = [] # 
        self.Kpp = []
        self.KtT = []
        self.L = []
        self.theta = [sigma, l, sigma_noise]
        self.theta0 = self.theta
        self.xmin = 0
        self.iter = 0
        self.padding = padding

        if len( x0 ) > 0:
            self.xmax = np.max(x0)
        else:
            self.xmax = 0
        Optimizer1D.__init__(self, tolerance, maxiter, procmax, x0, y0)
        
    def addpoints(self, newx, newy):
        
        if ( type(self.x) is np.ndarray):
            self.x = np.ndarray.tolist(np.squeeze( self.x ))
            self.y = np.ndarray.tolist(np.squeeze( self.y ))
            
        
        Optimizer1D.addpoints(self, newx, newy)
        
        self.x = np.reshape(self.x, (len(self.x), 1) )
        self.y = np.reshape(self.y, (len(self.y), 1) )
        
    def updateTrials(self):
        '''Update trial points depending on the set of training data.'''
        
        pp_interval = 10
        
        self.xt = np.array([])
        xs = np.sort(np.squeeze(self.x))
        for i in range(len(xs) - 1 ):
            self.xt = np.append(self.xt, np.linspace(xs[i],xs[i+1], pp_interval, endpoint = False) )
        self.Nx = len( np.squeeze(self.xt) )
        self.xt = np.reshape( self.xt, (self.Nx, 1) )

        self.Kpp = self.kernel(self.xt, self.xt, self.theta,measnoise = 0.)
        
    def logPosterior(self, theta, *args):
        '''
        log of the posterior probability as a function of the hyper
        parameters theta = [sigma, l, sigma_noise].
        '''
        x, y = args
        try:
            k = self.kernel( x, x, theta)
            L = np.linalg.cholesky(k)
        except( np.linalg.LinAlgError ):
            #print "WARNING: Linalg error, trying different theta..."
            for d in range(len(theta)):
                theta[d] = self.theta0[d]
                #print theta
        error = ""
        trymax = 100
        for i in range( trymax ):
            try:
                k = self.kernel( x, x, theta)
                L = np.linalg.cholesky(k)
                break
            except np.linalg.LinAlgError as e:
                error = str( e )
                #print "WARNING: Linalg error, trying different theta..."
                for d in range(len(theta)):
                    theta[d] += 0.5
                    #print theta
                    #print self.maxloc
        if i == trymax-1:
            print "FATAL, could not resolve linalg error: " + error
            print "x = " + str( np.sort(x.transpose()) )
            return 0
                    
        beta = np.linalg.solve(np.transpose(L), np.linalg.solve(L,y))
        logp = -0.5*np.dot( y.transpose(), beta ) - np.sum( np.log( np.diag(L) ) ) - np.shape(x)[0] / 2. * np.log( 2*np.pi )
        
        return -logp[0]

    def gradLogPosterior(self, theta, *args):
        '''
        gradient of the log of the posterior probability, with respect to the
        paratmeres theta = [sigma, l, sigma_noise].
        '''
        x, y = args
        d = len(theta)
        
        try:
            K = self.kernel(x,x,theta,wantderiv=True)
            L = np.linalg.cholesky(np.squeeze(K[:,:,0]))
        except( np.linalg.LinAlgError ):
            #print "WARNING: Linalg error, trying different theta..."
            for dd in range(len(theta)):
                theta[dd] = self.theta0[dd]*0.1
                #print theta
        error = ""
        trymax = 100
        for i in range(trymax):
            try:
                K = self.kernel(x,x,theta,wantderiv=True)
                L = np.linalg.cholesky(np.squeeze(K[:,:,0]))
                break
            except np.linalg.LinAlgError as e :
                error = str( e )
                #print "WARNING: Linalg error, trying different theta..." + str( e )
                for dd in range(len(theta)):
                    theta[dd] += 0.1*self.theta0[dd]
                    #print theta
                    #print self.maxloc
        if i == trymax-1:
            print "FATAL, could not resolve linalg error: " + error
            print "x = " + str( np.sort(x.transpose()) )
            return np.zeros(d)
            
            
        invk = np.linalg.solve(L.transpose(),np.linalg.solve(L,np.eye(np.shape(x)[0])))
    
        dlogpdtheta = np.zeros(d)
        for dd in range(1,len(theta)+1):
            
            dlogpdtheta[dd-1] = 0.5*np.dot(y.transpose(), np.dot(invk, np.dot(np.squeeze(K[:,:,dd]), np.dot(invk,y)))) - 0.5*np.trace(np.dot(invk,np.squeeze(K[:,:,dd])))
           

        #print dlogpdtheta
        return -dlogpdtheta
        
    def evalmeancov(self):
        '''Evaluates the mean and covariance of the GP.'''
        
        self.updateTrials()
        K = self.kernel(self.x, self.x, self.theta)
        
        try:
            L = np.linalg.cholesky(K)
            beta = np.linalg.solve(np.transpose(L), np.linalg.solve(L,self.y))
            #invk = np.linalg.solve(L.transpose(),np.linalg.solve(L,np.eye(np.shape(self.x)[0], dtype=np.double)))
            
        except np.linalg.linalg.LinAlgError as e:
            print "linalg error: " + str( e )
            return -1
        
        Kt = self.kernel(self.x, self.xt, self.theta, measnoise = 0.)
        
        Kt = np.array(Kt)
        Kpp = np.array(self.Kpp)
        KtT = np.transpose(Kt)
        
        self.mean = np.dot(KtT,beta)
        #self.mean = np.dot(KtT,np.dot(invk,self.y))
        
        v = np.linalg.solve(L, Kt)
        self.cov = Kpp - np.dot(np.transpose(v),v)
        #self.cov = Kpp - np.diag( np.dot( KtT, np.dot( invk, Kt ) ) )
        
        self.u = self.util( np.squeeze( self.mean ), self.cov )
        
        self.umax = np.max(self.u)
        maxloc = np.argmax(self.u)

        self.maxloc = self.xt[maxloc]
        
        
    def nextstep(self):
        
        args = (self.x, self.y)
        
        xs = np.sort(np.squeeze(self.x))
        
        xmin = np.max(xs)
        xmax = np.min(xs)
        xav = 0.
        for i in range(len(xs)-1):
            xtry = np.abs( xs[i]-xs[i+1] )
            xav += xtry
            if xtry < xmin and xtry > 1e-10:
                xmin = xtry
            if xtry > xmax:
                xmax = xtry
        xav /= float( len(xs) )
        
        #self.theta0[1] = (xmax + xmin)/2.
        self.theta0[1] = xav
                                  
        #print "xmin, xmax, xav = ",xmin,xmax,xav
        
        print "old theta = ", self.theta0, self.logPosterior(self.theta0,*args)
        newtheta = so.fmin_cg(self.logPosterior, self.theta0, fprime=self.gradLogPosterior, args=args, gtol=1e-4, maxiter=100, disp=1)
        print "new theta = ", newtheta, self.logPosterior(newtheta,*args)
        self.theta = newtheta
        
        # build the covariance matrix:
        
        self.evalmeancov()
        
        self.u = self.util(np.squeeze( self.mean ), self.cov)
        
        self.umax = []
        maxloc = []
       
        for iproc in range(0, self.pmax):
            uproc = []
            xproc = []
            while len(xproc) == 0:
                for it in range( len(self.xt) ):
                    xt = self.xt[it]
                    add = True
                    
                    for xp in np.append( self.x , maxloc):
                        if( np.abs( xt - xp ) < self.padding*self.tol and iproc > 0):
                            add = False
                            break
                    if add:
                        uproc.append( self.u[ it ] )
                        xproc.append( self.xt[it] )
                self.padding -=0.5
                if self.padding <= 0:
                    self.padding = 10
            xadd = xproc[np.argmax( uproc )]
            if xadd not in np.append( self.x, maxloc ):
                self.umax.append( np.max(uproc ))
                maxloc.append( xadd )
            else:
                if( iproc == 0):
                    print "Converged for x, ix, u(x) = ", xadd, np.argmax( uproc ), np.max( uproc )
                    self.converged = 1
                else:
                    exit()

        self.maxloc = maxloc
        
        if(len( self.maxloc ) == 0):
            return -1
        
        # Check convergence:
        for i in range(len(self.x)):
            for j in range(i+1, len(self.x)):
                if np.abs( self.x[i] - self.x[j] ) < self.tol:
                    self.converged = 1
                    print "Converged!"
                    print self.converged
                    break
            if self.converged:
                break
        
        return self.maxloc
    
    def minimize(self, model, sgenerator, pathwd, pathresults=None):
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
            model.gatherResults(sgenerator.structures[-len(newx):], pathwd, pathresults = pathresults, runprog = True)
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
    
    def minimize_parameters(self, model, hutil, plot = False):
        
        self.iter = 0
        while self.converged == 0:
            self.iter += 1
            
            if (self.iter > self.maxits):
                self.converged = -1
                return self.converged
            
            newx = self.nextstep()
            if self.converged:
                print self.x[np.argmin(self.y)]
                break
            print newx
            if(np.shape(newx) == np.shape(1) ):
                if newx == -1:
                    print "Update failed! Stopping."
                    self.converged = -1
                    break
            for xx in newx:
                if xx > hutil.imax-2:
                    xx = hutil.imax - 2
                elif xx < 2:
                    xx = 2
            newy = []
            for xx in newx:
                newy.append( -float(model.getMerit(hutil.interp_coords_from_dist(xx)/float(2**hutil.p))) )
        
            xstr = ''
            for xx in self.x:
                xstr += str( xx[0] ) + ', '
            #print "x = np.array([[ " + xstr + "]])"
            
            ystr = ''
            for yy in self.y:
                ystr += str( yy[0] ) + ', '
            #print "t = np.array([[ " + ystr + "]])"
        
            self.addpoints(newx,newy)
            if(plot):
                self.plotGP(model, hutil)
            
            
            
            
            print "Converged: ", self.converged
            #print "x = ",self.x
            #print "y = ",self.y

        return self.converged
        
    def plotGP(self, model, hutil):
        '''Plot the GP at its current stage in the minimization.
        
        Parameters
        
        model: The actual model- '''
        
        xt = self.xt; mean = np.squeeze( self.mean ); cov = self.cov; u = self.u;
        maxloc = self.maxloc; umax = self.umax
        x = self.x; y = self.y;
        
        var = np.reshape( np.abs( np.diag( self.cov ) ), (self.Nx,1) )
        
        yt = []
        [yt.append(-model.getMerit((hutil.interp_coords_from_dist(xx)/float(2**hutil.p)))) for xx in xt]
            
        pl.figure(5)
        pl.hold(False)
        pl.plot(xt,mean)
        pl.hold(True)
        pl.plot(x,y,'*')
        pl.fill_between(np.squeeze(xt), mean-np.squeeze(2*np.sqrt(var)),mean+np.squeeze(2*np.sqrt(var)), facecolor = "grey", alpha=0.5)
        pl.plot(xt,yt,'-')
        pl.plot(xt,u/5,'-g')
        pl.plot(maxloc,umax,'*')
        pl.ylim(-5, 5)
        
        lim_diff = (np.max(y) - np.min(y))/2.
        
        pl.gca().set_xlim( 0, np.max(x) )
        pl.gca().set_ylim( np.min(y) - lim_diff, np.max(y) + lim_diff )
        
        #writer.grab_frame()
        pl.pause(0.01)

    def getbest(self):
        index = np.argmin(self.y)
        return [self.x[index], self.y[index] ]
        
    def util(self, mean, cov):
            #return np.multiply(-mean + 0.001*np.diag(cov),100*np.diag(cov))
            #return -np.multiply(mean,np.diag(cov)) + 0.2*np.diag(cov)
            return np.exp(-mean) * ( 1 + np.sqrt( np.abs( np.diag(cov) ) ) )
            #return -mean
        
    def kernel(self,data1,data2,theta,wantderiv=False,measnoise=1.):
        '''Author: Stephen Marsland, 2014. Modified by Martin Franckie 2018'''
        
        theta = np.squeeze(theta)
        
        if (len(theta) == 2):
            theta = np.append(theta, 0.)
        
            
        # Gaussian
        if np.shape(data1)[0] == len(data1):
            d1 = np.shape(data1)[0]
            n = 1
        else:
            (d1,n) = np.shape(data1)
        d2 = np.shape(data2)[0]
        sumxy = np.zeros((d1,d2), dtype=np.double)
        for d in range(n):
            D1 = np.transpose( [data1[:,d]]) * np.ones((d1,d2), dtype=np.double )
            D2 = [data2[:,d]] * np.ones((d1,d2), dtype = np.double)
            sumxy += (D1-D2)
    
        k = theta[0]**2. * np.exp(-sumxy**2./(2.0*theta[1]**2.))
    
        if wantderiv:
            K = np.zeros((d1,d2,len(theta)+1), dtype=np.double)
            K[:,:,0] = k + measnoise*theta[2]**2.*np.eye(d1,d2, dtype=np.double)
            K[:,:,1] = 2.0 *k /theta[0]
            K[:,:,2] = k*sumxy**2./theta[1]**3.
            K[:,:,3] = 2.0*theta[2]*np.eye(d1,d2, dtype=np.double)
            return K
        else:    
            return k + measnoise*theta[2]**2.*np.eye(d1,d2, dtype = np.double)
        
    def kernel2(self,data1,data2,theta,wantderiv=False,measnoise=0):
        '''
        Author: Stephen Marsland, 2014 
        '''
        # Uses exp(theta) to ensure positive hyperparams
        theta = np.squeeze(theta)
        theta = np.exp(theta)
        # Squared exponential
        if np.ndim(data1) == 1:
            d1 = np.shape(data1)[0]
            n = 1
            data1 = data1*np.ones((d1,1))
            data2 = data2*np.ones((np.shape(data2)[0],1))
        else:
            (d1,n) = np.shape(data1)
    
        d2 = np.shape(data2)[0]
        sumxy = np.zeros((d1,d2))
        for d in range(n):
            D1 = np.transpose([data1[:,d]]) * np.ones((d1,d2))
            D2 = [data2[:,d]] * np.ones((d1,d2))
            sumxy += (D1-D2)**2*theta[d+1]
    
        k = theta[0] * np.exp(-0.5*sumxy) 
        #k = theta[0]**2 * np.exp(-sumxy/(2.0*theta[1]**2)) 
    
        #print k
        #print measnoise*theta[2]**2*np.eye(d1,d2)
        if wantderiv:
            K = np.zeros((d1,d2,len(theta)+1))
            K[:,:,0] = k
            K[:,:,1] = k
            K[:,:,2] = -0.5*k*sumxy
            #K[:,:,3] = theta[2]*np.eye(d1,d2)
            return K
        else:    
            return k
        
    def testopt(self):
        #theta = np.array([0.05,2]) # GP4
        #x = np.array([[-3.5, -2.5, -.5, .4, 2.25]]).transpose()
        #t = 0.55*np.array([[-2., 0., 1., 2., -1.]]).transpose()
        
        #x =  np.array([[    0., 2047.5,4095. , 3085.61122244, 2921.48296593, 1050.42084168,976.56312625, 3643.64729459, 3184.08817635, 1600.250501, 1378.67735471, 361.08216433, 287.2244489, \
        # 2461.9238477,         1854.6492986,         3947.28456914,         2896.86372745]]).transpose()
        #t =  np.array([[ 7.62307218, 2.53109568,          0.66060117,          0.7115907 ,          1.33369934,          1.41822756,          1.72123589,          2.16406453,\
        #  0.38066725,          1.19941393,          2.22988072,          0.90223427,          3.76357383,          1.26595732,          0.44283836,          3.1865146 ,\
        #  0.7630315 ]]).transpose()
        
        x = np.array([[ 0.0, 2047.5, 4095.0, 1012.47747748, 807.522522523, 3144.00900901, 2631.62162162, 1537.16216216, 3590.81081081, 3348.96396396, 2984.14414414, 2348.78378378, 327.927927928, 3197.2972973, 3836.75675676, 1278.91891892, 135.27027027, 2791.48648649, ]])
        t = np.array([[ 1.48640289174, 1.84806868376, 2.86214631467, 1.25035441867, 1.5632464401, 0.699012694924, 1.03543082705, 3.50768230106, 0.98214564022, 0.794669619712, 0.882179230782, 1.90127956833, 2.53058051441, 0.563755045005, 2.17131537249, 1.06258374857, 1.84409674801, 1.82407773281, ]])

        deltaxmin = np.max(x)
        deltaxmax = np.min(x)
        
        Nx = np.shape(x)[1]
        
        x = x.transpose()
        t = t.transpose()
        
        xs = np.sort(np.squeeze(x))
        
        print "sort = ", xs
        
        for i in range(Nx-1):
            if np.abs( xs[i] - xs[i+1] ) < deltaxmin:
                deltaxmin = np.abs(xs[i] - xs[i+1])
            if np.abs( xs[i] - xs[i+1] ) > deltaxmax:
                deltaxmax = np.abs(xs[i] - xs[i+1])
            
        
        
        print deltaxmin, deltaxmax
        
        #theta = np.array([4,2*deltaxmin]) # GP4
        theta = self.theta
        
        print "x = ", x
        Nx = 1000
        args = (x,t)
        print theta, -self.logPosterior(theta,*args)
        newTheta = so.fmin_cg(self.logPosterior, theta, fprime=self.gradLogPosterior, args=args, gtol=1e-4,maxiter=100,disp=1)
        print newTheta, -self.logPosterior(newTheta,*args)
        theta = newTheta
    
        xstar = np.reshape(np.linspace(np.min(x),np.max(x),Nx),(Nx,1))
    
        k = self.kernel(x,x,theta, measnoise = 1.)
        #kstar = [self.kernel(x,xs*np.ones((1,1)),self.theta,wantderiv=False) for xs in xstar]
        #kstar = np.squeeze(kstar)
        kstar = self.kernel(x,xstar, theta, measnoise = 0.)
        #kstarstar = [self.kernel(xs*np.ones((1,1)),xs*np.ones((1,1)),self.theta,wantderiv=False) for xs in xstar]
        #kstarstar = np.squeeze(kstarstar)
        #kstarstar = kernel2(xstar,xstar,theta,wantderiv=False)
        kstarstar = self.kernel(xstar, xstar, theta, measnoise = 0.)
    
        L = np.linalg.cholesky(k)
        invk = np.linalg.solve(L.transpose(),np.linalg.solve(L,np.eye(np.shape(x)[0])))
        #invL = np.linalg.inv(L)
        #invk = np.dot(invL.T,invL)
        print np.shape(kstar), np.shape(invk), np.shape(t)
        mean = np.dot(np.transpose(kstar),np.dot(invk,t))
        #print np.shape(kstarstar), np.shape(kstar), np.shape(invk)
        var = kstarstar - np.diag(np.dot(np.transpose(kstar),np.dot(invk,kstar)))
        print np.shape(mean), np.shape(var)
        
        #var = kstarstar - np.dot(kstar.transpose(),np.dot(invk,kstar))
        var = np.reshape( np.abs( np.diag( var ) ), (Nx,1) )
        print np.shape(var)
        #print mean
    
        pl.figure()
        pl.plot(xstar,mean,'-k')
        #pl.plot(xstar,mean+2*np.sqrt(var),'x-')
        #pl.plot(xstar,mean-2*np.sqrt(var),'x-')
        #print np.shape(xstar), np.shape(mean), np.shape(var)
        pl.fill_between(np.squeeze(xstar),np.squeeze(mean-2*np.sqrt(var)),
                        np.squeeze(mean+2*np.sqrt(var)),color='0.75')
        pl.plot(x,t,'ko')
        pl.axis('tight')
        pl.xlabel('x')
        pl.ylabel('f(x)')
        pl.show()
        
    def writeresults(self, pathresults, filename):
        with open(pathresults + "/" + filename, 'w') as f:
            [f.write( str( np.squeeze( self.x[i] )) + " " + str( 
                np.squeeze( self.y[i] )) +"\n") for i in range(0,len(self.x))]