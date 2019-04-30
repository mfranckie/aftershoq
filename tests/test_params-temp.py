from aftershoq.materials import *


gaas300 = GaAs(T=300)
gaas0 = GaAs(T=0)
print("GaAs:\n")
for key in gaas300.params:
    print(f"{key}: {gaas0.params[key]} (0K)")
    print(f"{key}: {gaas300.params[key]} (300K)")

alas300 = AlAs(T=300)
alas0 = AlAs(T=0)
print("\nAlAs:\n")
for key in alas300.params:
    print(f"{key}: {alas0.params[key]} (0K)")
    print(f"{key}: {alas300.params[key]} (300K)")

inas300 = InAs(T=300)
inas0 = InAs(T=0)
print("\nInAs:\n")
for key in inas300.params:
    print(f"{key}: {inas0.params[key]} (0K)")
    print(f"{key}: {inas300.params[key]} (300K)")

inp300 = InP(T=300)
inp0 = InP(T=0)
print("\nInP:\n")
for key in inp300.params:
    print(f"{key}: {inp0.params[key]} (0K)")
    print(f"{key}: {inp300.params[key]} (300K)")

GaSb300 = GaSb(T=300)
GaSb0 = GaSb(T=0)
print("\nGaSb:\n")
for key in GaSb300.params:
    print(f"{key}: {GaSb0.params[key]} (0K)")
    print(f"{key}: {GaSb300.params[key]} (300K)")

AlSb300 = AlSb(T=300)
AlSb0 = AlSb(T=0)
print("\nAlSb:\n")
for key in AlSb300.params:
    print(f"{key}: {AlSb0.params[key]} (0K)")
    print(f"{key}: {AlSb300.params[key]} (300K)")

InSb300 = InSb(T=300)
InSb0 = InSb(T=0)
print("\nInSb:\n")
for key in InSb300.params:
    print(f"{key}: {InSb0.params[key]} (0K)")
    print(f"{key}: {InSb300.params[key]} (300K)")

AlGaAs300 = AlGaAs(x = 0.25, T=300)
AlGaAs0 = AlGaAs(x = 0.25, T=0)
print("\nAlGaAs:\n")
for key in AlGaAs300.params:
    print(f"{key}: {AlGaAs0.params[key]} (0K)")
    print(f"{key}: {AlGaAs300.params[key]} (300K)")
