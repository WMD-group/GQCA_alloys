"""
This program calculates the binodal and spinodal
points of the phase diagram of a pseudobinary
alloy A_{1-x}B_{x}C

Clovis Caetano - July/2015
"""
import numpy as np
from scipy.optimize import brentq
from scipy.optimize import minimize_scalar
# The GQCA free energy is used, but any other expression also could be used
from gqca2 import free

# Variables
# x -> composition of the alloy
# T -> temperature (the unit depends on how the free energy is defined, but here is in kelvins)

# Cleaning the output file.
open('phase_diagram.dat', 'w').close()

# 2nd derivative of free energy (finite differences)
def d2free(x,T):
    h = 0.00001
    x1 = x-h
    x2 = x+h
    y = free(x,T)
    y1 = free(x1,T)
    y2 = free(x2,T)
    return (y2-2.0*y+y1)/h**2

# T0 -> Initial temperature
# dT -> Initial change of temperature
# The lowest temperature in the output will be T0 + dT
T0 = 150.0
dT = 50.0

# Initial guess to critical composition
xc = 0.5

# Itnitial guess to left spinodal point. It will be used to compute the new temperature.
xs0 = 0.0001

# Initial guesses for binodal points
x1_bin = 0.0001
x2_bin = 0.9999

# Headings of data to be shown
print "T (K)\tx1_bin\tx2_bin\tx1_spin\tx2_spin"

# Search will continue until the temperature variation be small enough
while (dT > 0.5):

    # Increasing the temperature
    T = T0 + dT

    # Calculating spinodal points using Brent method (d2free should be zero at these points)
    while True:
        try:
            x1_spin = brentq(d2free,x1_bin,xc,args=(T))
            x2_spin = brentq(d2free,xc,x2_bin,args=(T))
        # This is to prevent search above critical temperature
        except ValueError:
            # If it happens, dT is corrected
            T -= dT/2.0
            continue
        break

    # Critical composition is approximated by the average of spinodal points
    xc = (x1_spin+x2_spin)/2.0

    # Two auxiliary variables
    x1 = x1_spin
    x2 = x2_spin

    # These two functions are defined in order to find the common tangent
    def f1(x):
        return ((xc-x)*free(x2,T)+(x2-xc)*free(x,T))/(x2-x)
    def f2(x):
        return ((xc-x1)*free(x,T)+(x-xc)*free(x1,T))/(x-x1)

    # Search for binodal points
    h1 = 1.0
    h2 = 1.0
    k = 0
    while (h1>1e-6 or h2>1e-6):
        res = minimize_scalar(f1,bounds=(x1_bin,x1_spin), method='bounded')
        h1 = abs(x1-res.x)
        x1 = res.x
        res = minimize_scalar(f2,bounds=(x2_spin,x2_bin), method='bounded')
        h2 = abs(x2-res.x)
        x2 = res.x
        k += 1
        if k == 10:
            break
    # Binodal points
    x1_bin = x1
    x2_bin = x2

    # Showing the results
    print("%.1f\t%.4f\t%.4f\t%.4f\t%.4f" %(T,x1_bin,x2_bin,x1_spin,x2_spin))

    # Saving the results in output file
    with open("phase_diagram.dat","a") as output:
        output.write("%.2f %.5f %.5f %.5f %.5f\n" %(T,x1_bin,x2_bin,x1_spin,x2_spin))

    # Calculating the new temperature
    # The scale factor can be changed to increase the amound of points in the phase diagram
    dT = 0.01*(T-T0)/(x1_spin-xs0)
    T0 = T
    xs0 = x1_spin

# Printing the critical point
print "Critical composition = %.5f" %xc
print "Critical temperature = %.1f K" %T0
