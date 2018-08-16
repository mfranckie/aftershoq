'''
Created 2018

@author: Martin Franckie

Enum module for material parameters.
'''

meff = 0   # effective mass
Ec = 1     # Conduction band offset (eV)
Eg = 2     # Band gap (eV)
Ep = 3     # Kane energy (eV)
Valloy = 4 # Alloy scattering potential (eV)
ELO = 5    # Optical phonon energy (eV)
eps0 = 6   # Static dielectric constant
epsinf = 7 # High-frequency dielectric function
Vdef = 8   # Deformation potential (eV)
vlong = 9  # Longitudinal sound velocity (m/s)
massdens = 10   # Mass density (kg/m^3)
molV = 11  # Mol volume (m^3/mol)
lattconst = 12  # Lattice constant (A)
Nparam = 13     # Total number of parameters.

valdict = {"meff":meff, "CBO": Ec, "Eg": Eg, "Ep":Ep, "Alloy pot.": Valloy,
           "ELO": ELO, "eps(0)":eps0, "eps(inf)":epsinf,
           "deform. pot.":Vdef,"long. sound vel.":vlong,"mass dens":massdens,
           "mol volume":molV,"lattice constant": lattconst}

def initList(C):
    for _ in range(0,Nparam):
        C.append(0)