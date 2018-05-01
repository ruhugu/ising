#include "c_evolve.h"

int c_evolve_nofieldGlauber(
        int* spins_in, int* spins_out, int *neigh_list, int nspins,
        int* nneigh, int nneigh_max, double beta, long int nsteps)
{
    int j_MCstep;  // Counter for the MC step loop
    int j_step;  // Counter for the steps inside each MC step
    int j_neigh; // Counter for the neighbors loop
    int j_spin;
    int acceptprob_idx;
    double* acceptprob;
    long int naccept;

    // Initialize the number of accepted proposals to zero
    naccept = 0;

    // Allocate memory for the acceptance probability array
    acceptprob = (double*) malloc((2*nneigh_max + 1)*sizeof(double));

    // Precalculate the values of the Glauber acceptance probability
    // for the possible energy differences. Notice that the coupling
    // constant is set to 1).
    for (int j_boltz = 0; j_boltz <= 2*nneigh_max; j_boltz++)
    {
        double boltzfactor = exp(-beta*(2.*(j_boltz - nneigh_max)));
        acceptprob[j_boltz] = boltzfactor/(1. + boltzfactor);
    }

    // Store the initial value of the lattice in spins_out
    for (j_spin = 0; j_spin < nspins; j_spin++)
    {
        spins_out[j_spin] = spins_in[j_spin];
    }

    // MC steps Loop
    for (j_MCstep = 0; j_MCstep < nsteps; j_MCstep++)
    {
        // Each MC step consist in nspins spin updates, so 
        // that every spin "has the chance" to change.
        for (j_step = 0; j_step < nspins; j_step++)
        {
            // Choose a spin randomly
            j_spin = idran_(&nspins) - 1;

            // Calculate the index corresponding to the 
            // acceptance probability for this change.
            acceptprob_idx = 0;
            for (j_neigh = 0; j_neigh < nneigh[j_neigh]; j_neigh++)
            {
                acceptprob_idx += spins_out[
                        neigh_list[index2D(j_spin, j_neigh, nneigh_max)]];
            }
            acceptprob_idx = spins_out[j_spin]*acceptprob_idx + nneigh_max;
           
            // Accept the update with the corresponding
            // probability
            if (dranu_() < acceptprob[acceptprob_idx])
            {
                spins_out[j_spin] *= -1;
                naccept += 1;
            }
        }
    }

    // Free allocated memory
    free(acceptprob);
            
    return naccept;
}


int c_evolve_nofieldCouplGlauber(
        int* spins_in, int* spins_out, int *neigh_list, double* coupls,
        int nspins, int* nneigh, int nneigh_max, double beta, long int nsteps)
{
    int j_MCstep;  // Counter for the MC step loop
    int j_step;  // Counter for the steps inside each MC step
    int j_neigh; // Counter for the neighbors loop
    int j_spin;
    int neigh_idx;
    int aux_prob; // Auxiliary var to calculate acceptprob
    double acceptprob;
    double boltzfactor;
    long int naccept;

    // Initialize the number of accepted proposals to zero
    naccept = 0;

    // Store the initial value of the lattice in spins_out
    for (j_spin = 0; j_spin < nspins; j_spin++)
    {
        spins_out[j_spin] = spins_in[j_spin];
    }

    // MC steps Loop
    for (j_MCstep = 0; j_MCstep < nsteps; j_MCstep++)
    {
        // Each MC step consist in nspins spin updates, so 
        // that every spin "has the chance" to change.
        for (j_step = 0; j_step < nspins; j_step++)
        {
            // Choose a spin randomly
            j_spin = idran_(&nspins) - 1;

            // Calculate the the acceptance probability for this change.
            aux_prob = 0;
            for (j_neigh = 0; j_neigh < nneigh[j_neigh]; j_neigh++)
            {
                // Find the index in "spins_out" of the neighbour
                neigh_idx = neigh_list[index2D(j_spin, j_neigh, nneigh_max)];
                aux_prob += (
                        spins_out[neigh_idx]
                        *coupls[index2D(j_spin, neigh_idx, nspins)]);
            }
            aux_prob = spins_out[j_spin]*aux_prob;
            boltzfactor = exp(-beta*(2.*aux_prob));
            acceptprob = boltzfactor/(1. + boltzfactor);
           
            // Accept the update with the corresponding
            // probability
            if (dranu_() < acceptprob)
            {
                spins_out[j_spin] *= -1;
                naccept += 1;
            }
        }
    }
            
    return naccept;
}


double c_hamiltonian(int* spins, int* pairs, int npairs)
{
    double hamilt = 0;
    for (int j_pair=0; j_pair < npairs; j_pair++)
    {
        hamilt += - (
                spins[pairs[index2D(j_pair, 0, 2)]]
                *spins[pairs[index2D(j_pair, 1, 2)]]);
    }

    return hamilt;
}

// Given 2D indices, return the corresponding index in the flattened array
int index2D(int row, int column, int rowlength)
{
    return row*rowlength + column;
}            

// Find the maximum element in an int array
int max_element_int(int* array, int size)
{
    int maximum = array[0];

    for (int j = 1; j < size; j++)
    {
        if (array[j] > maximum)
        {
           maximum  = array[j];
        }
    }

    return maximum;
}

