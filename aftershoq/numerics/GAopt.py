import random
from aftershoq.numerics.optimizer import OptimizerGA
import numpy as np

class GAopt(OptimizerGA):

    def __init__(self, tolerance, maxiter, population, mutrate = 0.05,
                 mutsize = 0.1, corate = 1., nparents = None, noff = None,
                 limits = None,
                 x0 = [], y0 =  []):
        """Object for genetic algorithm optimization.

        Parameters
        ----------
        tolerance : float
            The tolerance for convergence. (TODO: explain how it is implemented)
        maxiter : int
            Maximum number of iterations until optimization is aborted.
        population : int
            Size of the population in each generation.
        mutrate : float
            Mutation rate. Gives the probability that a parameter is mutated.
        mutsize : float
            Mutation magnitude in percentage (valid for all parameters)
        corate : float
            Cross-over rate. Gives the probability that a parameter is used in
            the cross-over phase.
        nparents : int
            Number of parent pairs to be chosen.
        noff : int
            Number of offsprings for each pair of parents to produce.
        limits : array[(float, float)]
            Sets the parameter range limits, where each parameter has a limit
            defined by a tuple (min, max). If not set, no limits are imposed.
        x0 : array[ array [ float ] ]
            Initial paramter values.
        y0 : array[ float ]
            Intial fitness functions for each element of x0.
        """

        self.x, self.y = x0, y0
        self.tolerance  = tolerance
        self.maxiter = maxiter
        self.population = population
        self.mutrate = mutrate
        self.mutsize = mutsize
        self.corate = corate
        self.limits = limits
        if (nparents is None):
            self.nparents = 2
        else:
            self.nparents = nparents
        if (noff is None):
            self.noff = 2
        else:
            self.noff = noff


    def addpoints(self, newx, newy):
        """Add newly evaluated points.

        Parameters
        ----------
        newx : array[array[float]]
            Array of parameters for the newly evaluated points.
        newy : type
            Array of fitness values for the newly evaluated points.
        """

        for i in range(len(newx)):
            self.y.append(newy[i])
            self.x.append(newx[i])


    def nextstep(self):
        """Execute the evoloution for one generation:
        Step 1: Cross-over for current population.
        Step 2: Mutate new population.
        Step 3: Exclude the worst individuals from previous generation.

        Returns
        -------
        array[array[float]]
            Array of offspring parameter arrays, whose fitness have to be
            evalutated next.
        """

        newx = self.crossover()
        newx = self.mutate(newx)

        index_sort = np.argsort(self.y)

        oldy = self.y
        self.y = []
        oldx = self.x
        self.x = []
        Nsave = self.population - len(newx)
        for i in range(Nsave):
            self.y.append(oldy[ index_sort[i] ])
            self.x.append(oldx[ index_sort[i] ])

        return newx

    def mutate(self, params):
        """Executes the mutation phase of this iterations population.

        Returns
        -------
        array[array[float]]
            Array of new set of parameter values.

        """

        newx = []

        for xi in params:
            newxi = []
            rnd = random.random()
            for i in range(len(xi)):
                xij = xi[i]
                if self.limits is not None:
                    (lim_min, lim_max) = self.limits[i]
                if ( rnd < self.mutrate ):
                    # mutation happens
                    xmut = xij + self.mutsize*(random.random() - .5)*(lim_max - lim_min)
                    while self.limits is not None and \
                        (xmut < lim_min or xmut > lim_max):
                        D = lim_max - lim_min
                        xmut = xmut - np.floor((xmut - lim_min)/D)*D

                    newxi.append(xmut)
                else:
                    newxi.append(xij)

            newx.append(newxi)

        return newx


    def crossover(self):
        """Executes the crossover phase, in which parents in this generation
        are selected according to their fitness functions. The selected parents
        will generate one randomized offspring, then new parents are selected
        until the specified number of offsprings have been generated.

        Returns
        -------
        array[array[float]]
            Array of offspring parameter arrays.
        """

        index_sort = np.argsort(self.y)
        #print(f"Sorted fitness: {index_sort}")

        N = self.nparents
        Noff = self.noff

        #print(f"Crossover with {N} parents to {Noff} offspring.")

        #Z = (N*(N+3.) + 1.)/2.
        Z = 0.
        P = []
        for i in range( N ):
            Z += self.y[ index_sort[i] ]
            P.append( self.y[index_sort[i]] )
        P = np.array(P)/Z

        #for i in range(N):
        #    print(f"P_{i} = {P[i]}; m_{i} = {self.y[index_sort[i]]}")

        # Here we keep the already chosen parents, so that two are not chosen
        # twice.
        chosen_parents = np.zeros( (N,N) )

        newx = []



        for _ in range(self.nparents):

            p1, p2 = -1, -1

            # choose parents: p1 and p2 are indices of parents, sorted by decreasing
            # fitness function.
            while (p1 == -1 or p2 == -1):
                for i in range( N ):
                    rnd = random.random()*self.corate

                    #Pi = (N-i+1)/Z
                    Pi = P[i]


                    if (rnd < Pi):
                        if (p1 == -1):
                            p1 = i
                        else:
                            p2 = i

                    if (p2 != -1):
                        if chosen_parents[p1][p2] == 1:
                            continue
                        else:
                            break


            #print(f"Parents chosen: {p1} and {p2}.")
            chosen_parents[p1][p2] = 1

            # do the cross-over between p1 and p2:
            for _ in range(Noff):
                newxi = []
                x1 = self.x[ index_sort[p1] ]
                x2 = self.x[ index_sort[p2] ]
                for j in range( len( x1 ) ):
                    w = random.random() # each gene gets a random weight
                    co = x1[j]*w + x2[j]*(1-w)
                    newxi.append( co )

                newx.append(newxi)

        return newx
