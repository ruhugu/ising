import numpy as np
import random 
from matplotlib import pyplot as plt
from matplotlib import animation

import cevolve 


class Ising(object):
    """Base Ising model class.

    """
    def __init__(self, nspins, seed=None): 
        """Init method.

        Parameters
        ----------
            nspins : int
                Number of spins.

            seed : int
                Random number generator seed.

        """
        # Store parameters
        self.nspins = nspins

        # Create a flat array to store the spin configuration
        # (an structured array will be accesible through the
        # latt() method).
        self.spins = np.zeros(self.nspins, dtype="intc")

        # Set values randomly
        self.reset_random()

        # Create neighbour lists
        self.neighlist, self.nneigh = self.neighbours()

        # Initialize the random number generators
        if seed == None:
            seed = random.SystemRandom().randint(0, 32767)
        self.seed = seed  # Store seed
        cevolve.seed(seed)  # dranxor number generator
        np.random.seed(seed)  # Python random number generator
     

    # Methods to access the spins array with a structure.
    # These are placeholders to be replaced in the 
    # specific implementations of the Ising model.
    @property
    def latt(self):
        """Structured spin array.

        This is a placeholder.

        """
        return self.spins
    @latt.setter
    def latt(self, value):
        self.spins = value


    def neighbours(self):
        """Create neighbour list for the lattice.

        This is a placeholder.

        Returns
        -------
            neighbourlist : int array
                2D array with the neighbours of each cell.
                neighbourlist[i,j] is the index in self.latt
                of the j-th neighbour of the i-th cell (with
                i also its index in self.latt).

            n_neigh: int 
                Number of neighbours of each cell.
        
        """
        neighbourlist = None
        n_neigh = 0
        return neighbourlist, n_neigh


    def reset_random(self):
        """Set the value of the lattice spins randomly.

        """
        self.spins = np.random.choice([-1, 1], size=self.nspins).astype("intc")
        return 


    def evolve(self, nsteps, T):
        """Evolve the lattice in nsteps.

           Parameters
           ----------
               nsteps : int
                   Number of Monte Carlo steps.

               T : float
                   Temperature of the system in units of the inverse
                   of the Boltzmann's constant times the coupling
                   constant times the particle's spin squared.

        """
        # Raise an error if the given temperature is not valid
        if T <= 0:
            raise ValueError("The temperature must be greater than zero.")

        # Calculate the value of the thermodynamic beta
        # for the given temperature
        beta = 1./T

        # Call the time evolution C function
        self.spins = cevolve.evolve_nofieldGlauber(
                self.spins, self.neighlist, beta, nsteps)

        return


    # TODO: improve docs
    def slow_thermalization(self, T_final, T_ini=4., steps_per_T=1000,
                            n_T=30):
        """Thermalize the system slowly to the given temperature.

        Parameters
        ----------
            T_final : float
                Final temperature.

            T_ini : float
                Initial temperature in the thermalization process.
                
            steps_per_T : int
                Number of time steps per temperature.

            T_step : float
                Difference between consecutive temperature steps.

            n_T : int
                NUmber of temperatures.

        """
        # Create the temperature vector
        Ts = np.linspace(T_ini, T_final, n_T)

        # Temperature loop
        for T in Ts:
            self.evolve(steps_per_T, T)

        return


    def magnetization(self):
        """Calculate the instantaneous mean magnetization.

        """
        return np.mean(self.spins)


    def plot(self, size=3):
        """Plot the system configuration. 

        """

        fig, ax = plt.subplots(figsize=(size,size))
        im = ax.imshow(self.latt, cmap=self.cmap, vmin=-1, vmax=+1,
                       interpolation=None)
        return fig


    def animate(self, T, nframes=1000, steps_per_frame=1, frame_interval=300):
        """Animate the evolution of the lattice. 

        Parameters
        ----------
            T : float
                Temperature of the system.

            nframes : int
                Number of frames.

            steps_per_frame : int
                Number of timesteps between frames.

            frame_interval : float
                Time in milliseconds between frames.

        Returns
        -------
            anim

        """

        def update(i, T, steps_per_frame, im, self):
            self.evolve(steps_per_frame, T)
            im.set_array(self.latt)
            return im

        fig, ax = plt.subplots()
        im = ax.imshow(self.latt, cmap="Greys", vmin=-1, 
                    vmax=+1, interpolation=None)
        cbar = fig.colorbar(im, ax=ax)

        anim = animation.FuncAnimation(fig, update, frames=nframes, 
                                       blit=False, interval=frame_interval,
                                       fargs=(T, steps_per_frame, im, self))
        return anim

      
class Ising2D(Ising):
    """2D Ising model lattice with periodic boundary conditions.

    """

    def __init__(self, ncols, nrows, seed=None):
        # Store parameters
        self.ncols = ncols
        self.nrows = nrows

        self.nspins = self.ncols*self.nrows
        Ising.__init__(self, self.nspins, seed=seed)
        

    @property
    def latt(self):
        """Structured spin array.

        """
        return np.reshape(self.spins, self.shape())
#    @latt.setter
#    def latt(self, value):
#        np.reshape(self.spins, self.shape()) = value


    def shape(self):
        """Return lattice shape.

        """
        return [self.nrows, self.ncols]

    def neighbours(self):
        """Create neighbour list for the lattice.

        Each spin has its the closest ones (up, down, left, right)
        has neighbours, with periodic boundary condiions taked 
        into account.

        Returns
        -------
            neighlist : int array
                2D array with the neighbours of each cell.
                neighbourlist[i,j] is the index in self.latt
                of the j-th neighbour of the i-th cell (with
                i also its index in self.latt). The j values
                correspond, in increasing order, to the
                top, right, bottom and left neighbours.

            nneigh: int 
                Number of neighbours of each cell.
        
        """
        # Number of neighbours per spin.
        nneigh = 4 
        neighlist = np.zeros((self.nspins, nneigh), dtype="intc")

        for j_spin in range(self.nspins):
            # Calculate the row and column of the spin.
            [row, col] = np.unravel_index(j_spin, self.shape())

            # Store the neighbours' indexes.

            # Use the wrap parameter to implement the periodic boundaries.
            neighlist[j_spin][0] = np.ravel_multi_index(
                    [row-1, col], self.shape(), mode="wrap") # Top neighbour
            neighlist[j_spin][1] = np.ravel_multi_index(
                    [row, col+1], self.shape(), mode="wrap") # Right neighbour
            neighlist[j_spin][2] = np.ravel_multi_index(
                    [row+1, col], self.shape(), mode="wrap") # Bottom neighbour
            neighlist[j_spin][3] = np.ravel_multi_index(
                    [row, col-1], self.shape(), mode="wrap") # Left neighbour
            
        return neighlist, nneigh
