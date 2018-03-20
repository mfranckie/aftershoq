'''
Created on 15 Mar 2018

@author: martin
'''
from hilbert import HilbertCurve
import numpy as np

class HilbertUtil(object):
    '''
    Help class for interpolation between nodes of the hilbert curve
    '''


    def __init__(self, hilbert_curve):
        '''
        Constructor
        '''
        self.hilbert_curve = hilbert_curve
        self.N = hilbert_curve.n
        self.p = hilbert_curve.p
        self.imax = 2**(self.N*self.p)-1
        self.pmax = 2**self.p-1
        
    def interp_dist_from_coords(self,x):
        
        x = np.array(x)
        
        intx = np.floor(x)
        
        dx = x - intx
        
        # only one coordinate can move at a time = differ from an integer
        dd = np.max(dx)
        
        d = self.hilbert_curve.distance_from_coordinates(intx.astype(int)) + dd
        
        return d
    
    def interp_coords_from_dist(self,d):
        
        intd = min( int(d) , self.imax-1 )
            
        ddin = intd - d
        
        hmap1 = np.array( self.hilbert_curve.coordinates_from_distance( intd ) )
        hmap2 = np.array( self.hilbert_curve.coordinates_from_distance( intd+1 ) )
        
        dp = hmap1 - hmap2
        
        x = hmap1 + dp*ddin
        
        return x