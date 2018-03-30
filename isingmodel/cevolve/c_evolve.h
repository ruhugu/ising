#pragma once
#include "dranxor2/dranxor2C.h"
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* void c_evolve_nofieldGlauber
Evolves the system with in the given number of Monte Carlo steps.

The system is a 2D lattice with no external field and Glauber's 
method is used.

Parameters
----------
    spins_in : int pointer
        Initial state of the lattice. Must point to an array with 
        length nspins.

    spins_out : int pointer
        The final state of the lattice is saved here. Must point 
        to an array with length nspins.

    neigh_list : int pointer
        Pointer to a 2D array with the neighborsof each cell.
        neighborlist[i,j] is the index of the j-th neighbor of
        the i-th spin.

    nneigh : int
        Number of neighbors of each spin.

    nspins : int
        Number of spins in the lattice.

    beta : double
        Thermodynamic beta for the lattice temperature in units
        of the coupling constant.

    nsteps : long int
        Number of steps to evolve the system.

*/
void c_evolve_nofieldGlauber(
        int* spins_in, int* spins_out, int *neigh_list, int nspins,
        int nneigh, double beta, long int nsteps);


// Given 2D indices, return the corresponding index in the flattened array
int index2D(int row, int column, int rowlength);
