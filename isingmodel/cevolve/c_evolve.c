#include "c_evolve.h"

int c_evolve_nofieldGlauber(
        int* spins_in, int* spins_out, int *neigh_list, int nspins,
        int nneigh, double beta, long int nsteps)
{
    int j_MCstep;  // Counter for the MC step loop
    int j_step;  // Counter for the steps inside each MC step
    int j_neigh; // Counter for the neighbors loop
    int j_spin;
    int acceptprob_idx;
    long int naccept;
    double acceptprob[5];


    // Initialize the number of accepted proposals to zero
    naccept = 0;

    // Precalculate the values of the Glauber acceptance 
    // probability for the possible energy differences.
    // In this case, those are (setting the coupling
    // constant to 1): +8, +4, 0, -4, -8
    for (int j_boltz = 0; j_boltz<5; j_boltz++)
    {
        double boltzfactor = exp(-beta*(4*j_boltz - 8));
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
            for (j_neigh = 0; j_neigh < nneigh; j_neigh++)
            {
                acceptprob_idx += spins_out[
                        neigh_list[index2D(j_spin, j_neigh, nneigh)]];
            }
            acceptprob_idx = (spins_out[j_spin]*acceptprob_idx + 4)/2;
           
            // Accept the update with the corresponding
            // probability
            if (dranu_() < acceptprob[acceptprob_idx])
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

