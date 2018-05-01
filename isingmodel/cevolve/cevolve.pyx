# -*- coding: utf-8 -*-

# Python extension through Cython following:
# http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id10
# Nice example:
# https://gist.github.com/phaustin/4973792

from __future__ import (print_function, division, 
                        absolute_import, unicode_literals)

# Imports Cython declarations for numpy.
# cimport is used to import C data types, functions, etc defined in 
# other Cython file. Details: 
# http://cython.readthedocs.io/en/latest/src/userguide/sharing_declarations.html#the-cimport-statement
# In this case we are not importing Python numpy. We are loading the 
# Cython code that allows it to interact with numpy.
cimport numpy as cnp

# Enable Numpy-C-API access.
# Interesting:
# https://github.com/cython/cython/wiki/tutorials-numpy
# Numpy-C-API:
# https://docs.scipy.org/doc/numpy/reference/c-api.html
#np.import_array()
import numpy as np

# This tells Cython that there the following functions are defined 
# elsewhere and they header is in "c_evolve.h".
# cdef is used to define c functions.
# http://notes-on-cython.readthedocs.io/en/latest/function_declarations.html
# cdef extern especifies that the function is defined elsewhere.
# http://cython-docs2.readthedocs.io/en/latest/src/userguide/external_C_code.html
cdef extern from "c_evolve.h":
    int c_evolve_nofieldGlauber(
            int* spins_in, int* spins_out, int *neighlist, int nspins,
            int* nneigh, int nneigh_max, double beta, long int nsteps)
    int c_evolve_nofieldCouplGlauber(
            int* spins_in, int* spins_out, int *neigh_list, double* coupls,
            int nspins, int* nneigh, int nneigh_max, double beta,
            long int nsteps)
    double c_hamiltonian(int* spins, int* pairs, int npairs)

# More of the same for the random generator
cdef extern from "dranxor2/dranxor2C.h":
    void dranini_(int*)


# Define a wrapper function that will act as a bridge between Python and 
# the C function <--- no se hasta que punto es esto totalmente cierto
# not None: by default, Cython allows arguments that meet the specified
# data type or that are None. In order to prevent the last behaviour, we 
# must add not None after the parameter name.
# http://docs.cython.org/en/latest/src/userguide/extension_types.html#extension-types-and-none

# I can't find where the sintax np.ndarray[...] is explained. However, from 
# this example we can see that the first argument is the datatype, ndim 
# refers to the dimension of the array and I think mode determines how 
# the array is stored in the memory. In this case we would use the c-way
# (whatever that is). What I have found about the matter:
# Here they say it is realated to "efficient buffer access". Other modes are
# also presented.
# https://github.com/cython/cython/wiki/enhancements-buffer
# The <...> before argument are type casts:
# http://cython.readthedocs.io/en/latest/src/reference/language_basics.html#type-casting
# If we use a Python variable as an argument of a Cython function with a
# specified type, automatic conversion will be attempted.
# http://cython.readthedocs.io/en/latest/src/userguide/language_basics.html

# ALTERNATIVE: MemoryViews
# http://nealhughes.net/cython1/
# http://docs.cython.org/en/latest/src/userguide/memoryviews.html

def evolve_nofieldGlauber(
        cnp.ndarray[int, ndim=1, mode="c"] spins_in,
        cnp.ndarray[int, ndim=2, mode="c"] neighlist,
        cnp.ndarray[int, ndim=1, mode="c"] nneigh,
        double beta, int nsteps):

    nspins = spins_in.size

    # The C function mush know the number of columns in the array,
    # which coincides with the max number of neighbours
    nneigh_max = (neighlist.shape)[1]

    # Create array where the evolved configuration will be stored
    cdef cnp.ndarray[int, ndim=1, mode="c"] spins_out = \
            np.zeros(nspins, dtype="intc")

    naccept = c_evolve_nofieldGlauber(
            <int*> cnp.PyArray_DATA(spins_in),
            <int*> cnp.PyArray_DATA(spins_out),
            <int*> cnp.PyArray_DATA(neighlist.flatten()),
            nspins,
            <int*> cnp.PyArray_DATA(nneigh), nneigh_max, 
            beta, nsteps)

    return spins_out, naccept


def evolve_nofieldCouplGlauber(
        cnp.ndarray[int, ndim=1, mode="c"] spins_in,
        cnp.ndarray[int, ndim=2, mode="c"] neighlist,
        cnp.ndarray[double, ndim=2, mode="c"] coupls,
        cnp.ndarray[int, ndim=1, mode="c"] nneigh,
        double beta, int nsteps):

    nspins = spins_in.size

    # The C function mush know the number of columns in the array,
    # which coincides with the max number of neighbours
    nneigh_max = (neighlist.shape)[1]

    # Create array where the evolved configuration will be stored
    cdef cnp.ndarray[int, ndim=1, mode="c"] spins_out = \
            np.zeros(nspins, dtype="intc")

    naccept = c_evolve_nofieldCouplGlauber(
            <int*> cnp.PyArray_DATA(spins_in),
            <int*> cnp.PyArray_DATA(spins_out),
            <int*> cnp.PyArray_DATA(neighlist.flatten()),
            <double*> cnp.PyArray_DATA(coupls.flatten()),
            nspins,
            <int*> cnp.PyArray_DATA(nneigh), nneigh_max, 
            beta, nsteps)

    return spins_out, naccept


def hamiltonian(
        cnp.ndarray[int, ndim=1, mode="c"] spins,
        cnp.ndarray[int, ndim=2, mode="c"] pairs):

    # Calculate array length
    npairs = pairs.shape[0]

    hamilt = c_hamiltonian(
            <int*> cnp.PyArray_DATA(spins),
            <int*> cnp.PyArray_DATA(pairs.flatten()),
            <int> npairs)

    return hamilt


def seed(int iseed): 
    """Initialize the random number generator in the c extension.
    
    Parameters
    ----------
        iseed : int
            Seed for the random number generator.

    """
    dranini_(&iseed)
    return    
