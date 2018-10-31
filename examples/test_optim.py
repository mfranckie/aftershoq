'''
Created on 14 May 2018

@author: martin
'''
from aftershoq.interface import NDtestfunc
import numpy as np
from matplotlib import pyplot as pl
from aftershoq.numerics import Paraopt, Gaussopt
from aftershoq.utils import HilbertUtil
from aftershoq.hilbert_curve.hilbert import HilbertCurve
import random as rm

if __name__ == '__main__':
    
    print("Welcome to the test case for multi-variate optimization.")
    ND = 3
    NDconst = 0
    p = 6;
    
    loaddata = True
    
    # load calculated data:
    x0 = []
    y0 = []
    if loaddata:
        path_data = "/Users/martin/git/demo/results/hilbert.log"
        with open(path_data) as f:
            for line in f:
                x0.append(float( line.split()[0] )/1.)
                y0.append(float( line.split()[1] ))
    
    print(x0)
    print(y0)
    
    # the curse of dimensionality
    
    nvec = []
    volume = []
    volume.append( 2 )
    volume.append( np.pi)
    nvec.append(1)
    nvec.append(2)
    Nmax = 30
    for i in range( 3, Nmax ):
        volume.append( 2*np.pi/i*volume[i-3] )
        nvec.append(i)
        
    #pl.figure(2)
    #pl.plot(nvec, volume, '.-')
    #pl.xlabel('Dim')
    #pl.ylabel('Volume')
    #pl.show()
    
    scale = 2**p
    hc = HilbertCurve( p, ND+NDconst )
    hutil = HilbertUtil(hc)
    
    Ntest = 1
    Nits = 0
    Nconv = 0
    Nr = 1
    tol = 0.001
    Nproc = 10
    
    rmin = 0.3
    rmax = 2
    r_res = []
    rspace = np.linspace(rmin, rmax, Nr)
    
    print("Parameters:\nDim\tp\ttol\t#tests")
    print(str(ND) + "\t" + str(p) + "\t" + str(tol*hutil.imax) + "\t" + str(Ntest))
    print("imax = " + str( hutil.imax ))
    
    for r in rspace:
        Nits_r = 0
        Nconv_r = 0
        for itest in range(0, Ntest):
            print("test # ", itest)
            model = NDtestfunc(ND)
            
            #x0 = [0, (rm.random()*0.5+0.5)*hutil.imax, hutil.imax]
            #x0 = [0, 0.25*hutil.imax, 0.5*hutil.imax, 0.75*hutil.imax, hutil.imax]
        
            #y0 = []
        
            #[y0.append( -float(model.getMerit( hutil.interp_coords_from_dist(d)/float(scale) )) ) for d in x0]
            
            #x0 =  np.array([[    0., 2047.5,4095. , 3085.61122244, 2921.48296593, 1050.42084168,976.56312625, 3643.64729459, 3184.08817635, 1600.250501, 1378.67735471, 361.08216433, 287.2244489, \
            #                2461.9238477,         1854.6492986,         3947.28456914,         2896.86372745]]).transpose()
            #y0 =  np.array([[ 7.62307218, 2.53109568,          0.66060117,          0.7115907 ,          1.33369934,          1.41822756,          1.72123589,          2.16406453,\
            #                0.38066725,          1.19941393,          2.22988072,          0.90223427,          3.76357383,          1.26595732,          0.44283836,          3.1865146 ,\
            #                0.7630315 ]]).transpose()
            
            #print "Creating optimization obj."
            #opt = Paraopt(tolerance=tol*hutil.imax, maxiter=1000, procmax=Nproc, r=r, x0=x0, y0=y0)
            opt = Gaussopt(tolerance=tol*hutil.imax, maxiter=50, procmax=Nproc, 
                           x0=x0, y0=y0, 
                           sigma = 100, l = hutil.imax*0.01, sigma_noise=0.1,
                           padding=100, sigma_noise_max = 5)
            #opt = Gaussopt(tolerance=tol*hutil.imax, maxiter=50, procmax=Nproc, x0=x0, y0=y0, sigma = 2, l = 200, sigma_noise=0.1, padding=100)
            
            opt.evalmeancov()
            opt.plotGP()
            pl.show()
            opt.nextstep()
            #opt.evalmeancov()
            opt.plotGP()
            pl.show()
         

            #opt.testopt()
            #pl.show()
            #opt.minimize_parameters(model, hutil, plot=False)
            #opt.minimize_parameters(model, hutil)
            
            if(opt.converged and opt.iter > 3):
                Nits_r += opt.iter
                Nconv_r += 1
            
            print("Converged with code", opt.converged, "after", opt.iter, "iterations.")
        r_res.append(Nits_r/float(Ntest))
        Nits += Nits_r
        Nconv += Nconv_r
        
        #opt.plotGP(model, hutil)
        
    
    # ------ RESULTS --------
    
    print("Convergence rate: " + str( Nconv/float(Ntest*Nr) ))
    if Nconv == 0:
        print("No convergence reached.")
    else:
        print("Avg. iterations until convergence: " + str( Nits/float(Nconv) ))
    
    pl.figure(1)
    if( ND == 1):
        pl.plot(opt.x, opt.y)
        pl.show()
        exit(0)
    
    xres = []
    yres = []
    zres = []
    for newx in opt.x:
        par_res = hutil.interp_coords_from_dist(newx)/float(2**p)
        xres.append( par_res[0] )
        yres.append( par_res[1] )
        if ND > 2:
            zres.append( par_res[2] )
        
    
    pl.scatter(xres,yres, c = opt.y, marker='s', s=30)
    
    if ND > 2:
            pl.figure(2)
            print(np.shape(yres), np.shape(zres), np.shape(opt.y))
            pl.scatter(yres, zres, c = opt.y, marker='s', s=30)
    
    pl.figure(4)
    sortin = np.argsort(np.squeeze(opt.x), axis=0)
    xhil = []
    [xhil.append(opt.x[index]) for index in sortin]
    yhil = []
    [yhil.append(opt.y[index]) for index in sortin]
    pl.plot(xhil, yhil, '-*')
    pl.plot(x0,y0,'ro', hold=True)
    #pl.plot(opt.xt, opt.mean)
    
    
    
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
    if ND < 3:
        valpointz = np.append(valpointz, 0)
    print(valpointz)
    minval = 10000
    minvalz = 10000
    minz = x[2][0]
    miny = x[1][0]
    minx = x[0][0]
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
            
            if val[-1] < minval:
                minval = val[-1]
                minx   = x[0][i]
                miny   = x[1][j]
            if valz[-1] < minvalz:
                minvalz = valz[-1]
                miny   = x[1][j]
                minz   = x[2][i]
    
    mincoords = []
    [mincoords.append(x) for x in best]
    if ND < 3:
        mincoords.append(0)
    mincoords[0] = minx; mincoords[1] = miny; 
    mincoords[2] = minz
    mindist = hutil.interp_dist_from_coords(np.array(mincoords[0:ND])*float(scale))
    pl.plot(mindist, minval, '*', markersize=10, hold=True)
    pl.plot(mindist, minvalz, '*', markersize=10, hold=True)
    
    pl.figure(1)
    pl.scatter(xl,yl, c = val, marker = '+', s=40)
    pl.scatter(minx,miny,marker = '*', s=200, hold=True)
    
    best = hutil.interp_coords_from_dist(xbest)/float(scale)
    pl.scatter(best[0],best[1],marker='*', s=500, hold=True)
    pl.figure(2)
    pl.scatter(zl,yl, c = valz, marker = '+', s=40)
    pl.scatter(minz,miny,marker = '*', s=200, hold=True)
    minval = np.min(valz)
    bestval = np.min(val)
    print("Minimum = " + str( minval ) + ", optimization gives " + str( ybest ))
    
    best = hutil.interp_coords_from_dist(xbest)/float(scale)
    if ND < 3:
        best = np.append( best, 0 )
    pl.scatter(best[2],best[1],marker='*', s=500, hold=True)
    
    
    pl.figure(3)
    pl.plot(rspace,r_res)
    
    pl.show()
