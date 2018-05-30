'''
Created on 14 May 2018

@author: martin
'''
from interface.ndtestfunc import NDtestfunc
import numpy as np
from matplotlib import pyplot as pl
from numerics.paraopt import Paraopt
from utils.hilbertutil import HilbertUtil
from hilbert import HilbertCurve
import random as rm

if __name__ == '__main__':
    
    print "Welcome to the test case for multi-variate optimization."
    ND = 5
    NDconst = 0
    p = 6;
    
    
    scale = 2**p
    hc = HilbertCurve( p, ND+NDconst )
    hutil = HilbertUtil(hc)
    
    Ntest = 100
    Nits = 0
    Nconv = 0
    Nr = 10
    tol = 0.001
    Nproc = 10
    
    rmin = 0.3
    rmax = 2
    r_res = []
    rspace = np.linspace(rmin, rmax, Nr)
    
    print "Parameters:\nDim\tp\ttol\t#tests"
    print str(ND) + "\t" + str(p) + "\t" + str(tol) + "\t" + str(Ntest)
    
    for r in rspace:
        Nits_r = 0
        Nconv_r = 0
        for _ in range(0, Ntest):
            model = NDtestfunc(ND)
            
            x0 = [0, (rm.random()*0.5+0.5)*hutil.imax, hutil.imax]
        
            y0 = []
        
            [y0.append( -float(model.getMerit( hutil.interp_coords_from_dist(d)/float(scale) )) ) for d in x0]
            
            opt = Paraopt(tolerance=tol*hutil.imax, maxiter=1000, procmax=Nproc, r=r, x0=x0, y0=y0)
            
            opt.minimize_parameters(model, hutil)
            
            if(opt.converged):
                Nits_r += opt.iter
                Nconv_r += 1
        r_res.append(Nits_r/float(Ntest))
        Nits += Nits_r
        Nconv += Nconv_r
        
    
    # ------ RESULTS --------
    
    print "Convergence rate: " + str( Nconv/float(Ntest*Nr) )
    print "Avg. iterations until convergence: " + str( Nits/float(Nconv) )
    
    xres = []
    yres = []
    zres = []
    for newx in opt.x:
        par_res = hutil.interp_coords_from_dist(newx)/float(2**p)
        xres.append( par_res[0] )
        yres.append( par_res[1] )
        zres.append( par_res[2] )
    
    pl.scatter(xres,yres, c = opt.y, marker='s', s=30)

    pl.figure(2)
    pl.scatter(yres, zres, c = opt.y, marker='s', s=30)
    
    pl.figure(4)
    pl.plot(opt.x, opt.y)
    pl.plot(x0,y0,'ro', hold=True)
    
    
    [xbest, ybest] = opt.getbest()
    
    best = hutil.interp_coords_from_dist(xbest)/float(scale)
    best2 = hutil.interp_coords_from_dist(xbest)/float(scale)
    
    pl.plot(xbest,ybest,'*', markersize=20, hold=True)
    
    xmin = 0; xmax = 1
    x = []
    [ x.append( np.linspace(xmin, xmax, 50) ) for _ in range ( 0, 3 ) ]
    
    res = []
    xl = []
    yl = []
    zl = []
    val = []
    valz = []
    valpoint = best
    valpointz = best2
    for i in range (0, len(x[0])):
        res.append([])
        valpoint[0] = x[0][i]
        valpointz[2] = x[2][i]
        for j in range (0, len(x[1])):
            yl.append(x[1][j])
            xl.append(x[0][i])
            zl.append(x[2][i])
            
            valpoint[1] = x[1][j]
            valpointz[1] = x[1][j]
            res[i].append( -model.testfunc( valpoint ) )
            val.append( -model.testfunc( valpoint ) )
            valz.append( -model.testfunc( valpointz ) )
    
    pl.figure(1)
    pl.scatter(xl,yl, c = val, marker = '+', s=40)
    
    best = hutil.interp_coords_from_dist(xbest)/float(scale)
    pl.scatter(best[0],best[1],marker='*', s=500, hold=True)
    pl.figure(2)
    pl.scatter(zl,yl, c = valz, marker = '+', s=40)
    
    best = hutil.interp_coords_from_dist(xbest)/float(scale)
    pl.scatter(best[2],best[1],marker='*', s=500, hold=True)
    
    
    pl.figure(3)
    pl.plot(rspace,r_res)
    
    
    pl.show()