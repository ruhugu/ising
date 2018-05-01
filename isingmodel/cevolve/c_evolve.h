#pragma once
#include "dranxor2/dranxor2C.h"
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* int c_evolve_nofieldGlauber
Evolves the system with in the given number of Monte Carlo steps.

The system is a 2D lattice with no external field and Glauber's 
method is used. The couplings are assumed constant.

Parameters
----------
    spins_in : int array pointer
        Initial state of the lattice. Must point to an array with 
        length nspins.

    spins_out : int pointer
        The final state of the lattice is saved here. Must point 
        to an array with length nspins.

    neigh_list : int pointer
        Pointer to a 2D array with the neighborsof each cell.
        neighborlist[i,j] is the index of the j-th neighbor of
        the i-th spin.

    nneigh : int array
        Number of neighbors of each spin.

    nspins : int
        Number of spins in the lattice.

    beta : double
        Thermodynamic beta for the lattice temperature in units
        of the coupling constant.

    nsteps : long int
        Number of steps to evolve the system.

Returns
-------
    naccept : long int
        Number of accepted proposals.

*/
int c_evolve_nofieldGlauber(
        int* spins_in, int* spins_out, int *neigh_list, int nspins,
        int* nneigh, int nneigh_max, double beta, long int nsteps);


/* Evolves the system with in the given number of Monte Carlo steps.

The system is a 2D lattice with no external field and Glauber's 
method is used. Arbitrary couplings are allowed.

Parameters
----------
    spins_in : int array pointer
        Initial state of the lattice. Must point to an array with 
        length nspins.

    spins_out : int pointer
        The final state of the lattice is saved here. Must point 
        to an array with length nspins.

    neigh_list : int pointer
        Pointer to a 2D array with the neighborsof each cell.
        neighborlist[i,j] is the index of the j-th neighbor of
        the i-th spin.

    nneigh : int array
        Number of neighbors of each spin.

    nspins : int
        Number of spins in the lattice.

    beta : double
        Thermodynamic beta for the lattice temperature in units
        of the coupling constant.

    nsteps : long int
        Number of steps to evolve the system.

Returns
-------
    naccept : long int
        Number of accepted proposals.

*/
int c_evolve_nofieldCouplGlauber(
        int* spins_in, int* spins_out, int *neigh_list, double* coupls,
        int nspins, int* nneigh, int nneigh_max, double beta, long int nsteps);


// TODO: check how C handles the pointers to arrays, the docs may 
// be wrong

/* Evaluate the Hamiltonian of the system

 Parameters
 ----------
     spins : int array pointer  
         Spins in the lattice.

     pairs : int 2d array pointer
         2D array where pairs[i][j] is the j-th element (j=0,1) of
         the i-th pair.

     nparis : int
         Number of pairs.
         
  Returns
  -------
      hamilt : double
         Value of the Hamiltonian.

*/

double c_hamiltonian(int* spins, int* pairs, int npairs);

// Given 2D indices, return the corresponding index in the flattened array
int index2D(int row, int column, int rowlength);


// Find the maximum element in an int array
int max_element_int(int* array, int size);

