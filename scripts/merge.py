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
# File lists to be merged. The files will be merged element wise
filelist1 = ["isingR10C10.dat",
            "isingR20C20.dat",
            "isingR40C40.dat",
            "isingR80C80.dat",
            "isingR160C160.dat",
            ]

filelist2 = [#"isingR10C10_crit.dat",
            "isingR20C20_crit.dat",
            "isingR40C40_crit.dat",
            "isingR80C80_crit.dat",
            "isingR160C160_crit.dat",
            ]

# Filenames of the merged results
filelist_m = [#"isingR10C10_m.dat",
            "isingR20C20.dat",
            "isingR40C40.dat",
            "isingR80C80.dat",
            "isingR160C160.dat",
            ]


# Import data
results_list1 = list()
results_list2 = list()
for fname in filelist1:
    results_list1.append(isingmodel.Results(fname=fname))
for fname in filelist2:
    results_list2.append(isingmodel.Results(fname=fname))

# Merge
results_list_m = list()
for res1, res2 in zip(results_list1, results_list2):
    results_list_m.append(isingmodel.mergeresults([res1, res2]))

# Save to file 
for fname, res_m in zip(filelist_m, results_list_m):
    res_m.savetxt(fname=fname)
