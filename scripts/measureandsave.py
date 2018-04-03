# -*- coding: utf-8 -*-
# Measure and store lattices with different sizes and temperatures

import os
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import isingmodel


# Parameters
Ls = [10, 20, 40, 80, 160]  # Lattice lengths
#Ts = np.linspace(4, 0.1, 30)  # Temperatures
Ts = np.linspace(2.3, 2.15, 30)  # Temperatures
measureinterval = 1  # Number of steps between measures
relaxsteps = 1000  # Number of steps waited before measuring
nmeasures = 100000*np.ones(Ts.size, dtype=int)  # Number of measures per temp.



# Lattice size loop
for L in Ls:
    print("L: {0}".format(L))
    # Create the lattice and results object
    latt = isingmodel.Ising2D(L, L)
    results = isingmodel.Results((L,L))

    # Relax to the first temperature to be measured
    latt.slow_thermalization(Ts[0], T_ini=4., T_step=0.1, steps_per_T=1000)

    # Temperature loop
    for T_idx, (T, nmeasure) in enumerate(zip(Ts, nmeasures)): 
        # Let the lattice thermalize
        latt.evolve(relaxsteps, T)

        # Measure the system
        results.measure(T, nmeasure, latt, measureinterval=measureinterval)

    # Save data to file
    results.savetxt(fname="isingR{0}C{0}_crit.dat".format(L))
