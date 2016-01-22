"""
This program uses generalised quasichemical approximation (GQCA)
to calculate some thermodynamics properties of a pseudobinary
alloy A_{1-x}B_{x}C

Clovis Caetano -  July/2015
"""
import numpy as np
from scipy.optimize import brentq

# Boltzmann constant in eV per K (CODATA 2006)
kB = 8.617343e-5

# Variables
# x -> composition of the alloy
# T -> temperature

# Reading input data from file energies.dat
# j -> index of configurations
# nj -> number of atoms B in configuration j
# gj -> degeneracy of configuration j
# ej -> energy of configuration j in eV per formula unit
nj = []
gj = []
ej = []
with open("energies.dat","r") as input:
    for line in input:
        if line.startswith('#'):
            continue
        data = line.split()
        nj.append(int(data[1]))
        gj.append(int(data[2]))
        ej.append(float(data[4])/3) # Dividing by 3 in order to obtain energy/anion

# J -> total number of configurations
# N -> number of atoms in the sublattice
J = len(nj)
N = nj[J-1]

# Calculating the excess energies dej
dej = []
for j in range(J):
    x = float(nj[j])/float(N)
    dej.append(ej[j]-(1-x)*ej[0]-x*ej[J-1])

# This quantity is used in the definition of the poliynomium
def alpha(j,x,T):
    return (N*x-nj[j])*gj[j]*np.exp(-dej[j]/(kB*T))

# Poliynomium related to a constraint of GQCA probability
def polynom(eta,x,T):
    total = 0
    for j in range(J):
        total += alpha(j,x,T)*eta**nj[j]
    return total

# Brent method is used to find the root of the polynomial equation
def root(x,T):
    eta_0 = 0.0
    eta_1 = 1e0
    while polynom(eta_1,x,T) >= 0:
        eta_1 = 10*eta_1
    return brentq(polynom,eta_0,eta_1,args=(x,T),maxiter=1000)

# GQCA probability of finding the configuration j in the alloy
def xj(j,T,r0):
    total = 0
    for i in range(J):
        total += gj[i]*r0**nj[i]*np.exp(-dej[i]/(kB*T))
    return gj[j]*r0**nj[j]*np.exp(-dej[j]/(kB*T))/total

# Probability in a random alloy (a priori)
def xj0(j,x):
    return gj[j]*x**nj[j]*(1-x)**(N-nj[j])

# Enthalpy of mixing
def enthalpy(x,T):
    if x == 0 or x == 1:
        return 0
    else:
        r0 = root(x,T)
        total = 0
        for j in range(J):
            total += xj(j,T,r0)*ej[j]
        return total-(1-x)*ej[0]-x*ej[J-1]

# Enthalpy of mixing (random alloy)
def enthalpy0(x):
    if x == 0 or x == 1:
        return 0
    else:
        total = 0
        for j in range(J):
            total += xj0(j,x)*ej[j]
        return total-(1-x)*ej[0]-x*ej[J-1]


# Entropy of mixing
def entropy(x,T):
    if x == 0 or x == 1:
        return 0
    else:
        r0 = root(x,T)
        total = 0
        for j in range(J):
            total += xj(j,T,r0)*np.log(xj(j,T,r0)/xj0(j,x))
        return -kB*(x*np.log(x)+(1-x)*np.log(1-x)+total/N)

# Entropy of mixing (ideal solution)
def entropy0(x):
    if x == 0 or x == 1:
        return 0
    else:
        return -kB*(x*np.log(x)+(1-x)*np.log(1-x))

# Free energy of mixing (Helmholtz)
def free(x,T):
    return enthalpy(x,T)-T*entropy(x,T)
