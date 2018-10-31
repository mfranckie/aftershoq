'''
Created on 15 Mar 2018

@author: Martin Franckie

'''
from hilbert_curve.hilbert import HilbertCurve
import numpy as np

class HilbertUtil(object):
    '''
    Help class for interpolation between nodes of a HilbertCurve.
    '''


    def __init__(self, hilbert_curve):
        '''Constructor.
        hilbert_curve: A HilbertCurve instance.
        '''
        self.hilbert_curve = hilbert_curve
        self.N = hilbert_curve.n
        self.p = hilbert_curve.p
        self.imax = 2**(self.N*self.p)-1
        self.pmax = 2**self.p-1
        
    def interp_dist_from_coords(self,x):
        '''Interpolates and returns a distance along the curve from (unscaled)
        coordinates. The coordinates must lie on the hilbert curve.
        '''
        
        x = np.array(x)
        
        intx = np.floor(x)
        
        dx = x - intx
        
        # only one coordinate can move at a time = differ from an integer
        dd = np.max(dx)
        
        try:
            d = self.hilbert_curve.distance_from_coordinates(intx.astype(int)) + dd
        except( ValueError ):
            for i in range(len(intx)):
                while intx[i] > 2**self.p-1:
                    intx[i]-=1
            d = self.hilbert_curve.distance_from_coordinates(intx.astype(int)) + dd
        
        return d
    
    def interp_coords_from_dist(self,d):
        '''Interpolates and returns a set of (unscaled) coordinates from the
        given distance along the curve.
        '''
        
        intd = min( int(d) , self.imax-1 )
            
        ddin = intd - d
        
        hmap1 = np.array( self.hilbert_curve.coordinates_from_distance( intd ) )
        hmap2 = np.array( self.hilbert_curve.coordinates_from_distance( intd+1 ) )
        
        dp = hmap1 - hmap2
        
        x = hmap1 + dp*ddin
        
        return x
    
    def scale_coordinates(self, x, ranges):
        '''
        Returns the scaled coordinates on the interval [xmin,xmax], 
        where "ranges" contains the [min,max] values for each parameter in x.
        '''
        
        xscaled = []
        if np.ndim(x) == 1:
            d1 = np.shape(x)[0]
            d2 = 1
        else:
            (d1,d2) = np.shape(x)
        x = np.reshape(x, (d1, d2))
        for i in range(d1):
            xi = []
            for iparam in range(d2):
                xmin = ranges[iparam][0]
                xmax= ranges[iparam][1]
                xs = x[i][iparam]*(xmax-xmin)/float(self.pmax) + xmin
                xi.append(xs)
            xscaled.append(xi)
        return xscaled
    
    def unscale_coordinates(self, x, ranges):
        '''
        Returns the unscaled coordinates on the interval [0,pmax] , where
        "ranges" contains the [min,max] values for each parameter in x.
        '''
        xscaled = []
        if np.ndim(x) == 1:
            d1 = np.shape(x)[0]
            d2 = 1
        else:
            (d1,d2) = np.shape(x)
        x = np.reshape(x, (d1, d2))
        for i in range(d1):
            xi = []
            for iparam in range(d2):
                xmin = ranges[iparam][0]
                xmax= ranges[iparam][1]
                xs = (float(x[i][iparam])-float(xmin))/(float(xmax)-float(xmin))*float(self.pmax)
                xi.append(xs)
            xscaled.append(xi)
        return xscaled
