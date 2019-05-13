# test threading in Python

from aftershoq.qcls import EV2416
from aftershoq.structure import Sgenerator
from aftershoq.interface import Inegf
import threading
import time
import os

if __name__ == "__main__":

    s0 = EV2416()
    print(s0)
    sg = Sgenerator(s0,dw=[1,1,1,1,1,1])

    sg.genRanStructs(3)

    model = Inegf(binpath="/Users/martinfranckie/git/NEGFT8/bin/")

    model.merit = model.merits["max gain"]

    model.numpar["Nstates"] = 3
    model.numpar["NE"],model.numpar["Nk"] = 50,50
    model.numpar["maxits"] = 5
    model.numpar["Nefield"] = 5
    model.numpar["efield0"] = 0.050
    model.numpar["defield"] = 0.001
    print(model.numpar)

    path = os.getcwd() + "/sequence/"
    print("Running in " + path)

    #t = model.runSequence(sg.structures[0], path)

    numparIV = model.numpar.copy()
    numparGain = model.numpar.copy()
    numparGain["boolEins"] = True
    numparGain["Nh"] = 1
    numparGain["Nefield"] = 1
    numparGain["efac0"] = 1e-5
    numparGain["omega0"], numparGain["domega"], numparGain["Nomega"] = 0.010,0.001,5
    seq = [numparIV, numparGain]

    t = model.runStructSeq(sg.structures,path, seq = seq, runwannier=False)

    print("Done")

    for tt in t:
        print(f"Joining thread {t}")
        tt.join()

    model.gatherResults(sg.structures, path)

    for s in sg.structures:
        negft = model.getresults(structure=s,path=path,datpath="Gain")
        print(f"Results for sid={s.sid}: \n{negft}\nmeit={s.merit}")
