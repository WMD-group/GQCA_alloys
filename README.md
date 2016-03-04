# GQCA Alloy Analysis
A Python tool to analyse the thermodynamics of alloys through the Generalized Quasi-Chemical Approximation (GQCA) written by Clovis Caetano.

# Description
The GQCA approach to calculate the configurational free energy of alloys was first described in [A. Sher, M. van Schilfgaard, A. B. Chen and W. Chen, Physical Review B 36, 4279 (1987)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.36.4279) and well explained in [A. B. Chen and A. Sher, Semiconductor Alloys (Plenum Press, New York, 1995)](http://www.springer.com/us/book/9780306450525). The repository contains a Python script to perform analysis of thermodynamics of a solid solution as a statistical ensemble of independent structures. The analysis runs on the energy of the structures that describe the configurational space restricted to a specific crystallographic supercell.

*As input, the scripts requires:*

1) number of structures

2) number of species substituted with respect to pure end members

3) degeneracy of each structure

4) energy of the single cell

The code `gqca2.py` reads the file `energies.dat` and returns the enthalpy of mixing, configurational entropy of mixing and the Helmholtz free energy as functions of the alloy composition and temperature.

The script `phase_diagram.py` calls `gqca2.py` and calculates the binodal and spinodal points of the system for a given temperature. The code is made to automatically change the temperature, collect the phase diagram points and write them in the file `phase_diagram.dat`. If it doesn't work, it would be interesting to plot the free energy to understand its variation.

The script `plot_phase_diagram.py` just reads the file phase_diagram.py and draws the phase diagram in a pdf file. 
