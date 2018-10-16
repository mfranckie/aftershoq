'''
Created on 30 Jan 2018

@author: martin
'''
import time

class Numerics():
    '''
    classdocs
    '''


    def __init__(self, model):
        '''
        Constructor
        '''
        self.model = model
        self.processes = []
    
    def generateAndRun(self,gen,Ns,path):
        gen.genRanStructs(Ns)
        self.runStructures(gen.structures, path)
            
    def runStructures(self,structures,path):
        proc = self.model.runStructures(structures,path)
        [self.processes.append(i) for i in proc]
        
    def waitforproc(self,delay):
        pactive = True
        while pactive:
            print("Processes running...")
            pactive = False
            for p in self.processes:
                if self.model.pltfm.jobstatus(p):
                    pactive=True
                    #break
            time.sleep(delay)
        print("All processes terminated!")
