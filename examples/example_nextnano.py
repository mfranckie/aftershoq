from aftershoq.interface import Inextnano
from aftershoq.materials import GaAs, AlGaAs
from aftershoq.numerics.runplatf import Local
from aftershoq.qcls import EV2416
import aftershoq.utils.systemutil as su
import aftershoq.utils.debug as dbg

pathwd = "../demo/NextNano"
su.mkdir(pathwd)

dbg.open(dbg.verb_modes["chatty"], outfile=pathwd + "/debug.log")

nnroot = "C:/Users/Martin Franckie/Documents/nextnanoQCL_2018_03_29/nextnano/2018_03_29/nextnano.QCL/"

lic = nnroot + "../License/License_nnQCL.lic"

model = Inextnano(nnroot, Local(), GaAs(), lic)

s = EV2416()

for l in s.layers:
    if l.material.name == "AlGaAs":
        l.material.name = "Al(x)Ga(1-x)As"
print(s)

model.numpar["efield0"] = 0.050
model.numpar["defield"] = 0.005
model.numpar["Nefield"] = 5
model.numpar["NE"] = 400
model.numpar["Nk"] = 400
model.numpar["Nz"] = 100
# Override default (automatic calculation)
#model.numpar["Emax"] = 0.060
#model.numpar["Ekmax"] = 0.060

model.writeInputFile(s, pathwd)

exit(0)

proc = model.runStructures([s], pathwd)

[out, err] = proc[0].communicate()
print("output from program: \n" + out)
print("ERROR: \n" + err)

model.waitforproc(5, "nextnano.QCL is running...")
