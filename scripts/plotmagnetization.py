# -*- coding: utf-8 -*-

# Plot the magnetization versus the temperature with the data obtained
# with measureandsave.py.

import os
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import isingmodel


# Parameters
# Files where the data is stored
filelist = ["isingR10C10.dat",
            "isingR20C20.dat",
            "isingR40C40.dat",
            "isingR80C80.dat",
             "isingR160C160.dat"]
figylen = 3.5
figxlen = 1.8*figylen


# Import data
results_list = list()
for fname in filelist:
    results_list.append(isingmodel.Results(fname=fname))

# Find minimum and maximum temperatures in order to find the exact
# solution plot range
T_max = np.amax([np.amax(results.Ts) for results in results_list])
T_min = np.amin([np.amin(results.Ts) for results in results_list])


# Plot
fig, ax = plt.subplots(figsize=(figxlen, figylen))
ax.set_xlabel(r"$T$")
ax.set_ylabel(r"$m$")

# Data plot
for results in results_list:
    ax.errorbar(
            results.Ts, results.mags, yerr=results.mag_err(), fmt="o",
            capsize=2, markersize=4, alpha=0.7, 
            label="Shape={0}".format(results.shape))

# Onsager solution plot
Ts = np.linspace(T_min, T_max, 200)
ax.plot(Ts, isingmodel.Ising2D.mag_exact(Ts), "-",
            label="Exact solution", color="gray")

ax.legend()
fig.tight_layout()

# Show plot
plt.show()

# Save plot to file
#fig.savefig("magnetization.png", dpi=200)
