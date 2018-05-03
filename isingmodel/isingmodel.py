import sys
import os 

import numpy as np
import random 
from matplotlib import pyplot as plt
from matplotlib import animation

sys.path.insert(0, os.path.abspath(
        os.path.join(os.path.dirname(__file__), '../networks')))
import networks

import cevolve 


class Ising(object):
    """Base Ising model class.

    """
    def __init__(self, nspins, network=None, shape=None, seed=None): 
        """Init method.

        Parameters
        ----------
            nspins : int
                Number of spins.

            network : Network object
                Network object with the connection between spins.

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

        # If a network is not given create an empty one
        if network == None:
            self.network = networks.Network(self.nspins)
        else:
            self.network = network

        # If the lattice has no shape attribute, create it
        if shape == None:
            self._shape = tuple((self.nspins,))
        else: 
            if np.prod(shape) == self.nspins:
                self._shape = tuple(shape)
            else:
                raise ValueError(
                        "The given shape does not match the number of spins.")
            
        # Create neighbour lists
        self.update_neighbours()

        # Initialize the random number generators
        if seed == None:
            seed = random.SystemRandom().randint(0, 32767)
        self.seed = seed  # Store seed
        cevolve.seed(seed)  # dranxor number generator
        np.random.seed(seed)  # Python random number generator


    # Methods to access the spins array with a structure.
    @property
    def latt(self):
        """Structured spin array.

        """
        return np.reshape(self.spins, self.shape())
    # TODO: create a setter method that works with shape
#    @latt.setter
#    def latt(self, value):
#        self.spins = value

    def shape(self):
        """Return lattice shape.

        """
        return self._shape


    def update_neighbours(self):
        """Update the lattice structure list using network attribute.

        """
        # Calculate the number of neighbours of each spin
        ns_neigh = self.network.degree_out.astype("intc")
        nneigh_max = np.amax(ns_neigh)

        self.nneigh = ns_neigh

        # Store the pairs of neighbouring spins
        # (this is useful for calculating the hamiltonian)
        self.neighpairs = np.array(
                self.network.adjlist(directed=False, weighted=False),
                dtype="intc")
        
        # Initialize neighbour index and number lists.
        # If a spin has less than nneigh_max neighbours, the first columns
        # of its row are filled with its neighbours and the rest with -1's
        # TODO: is there are spins with much more neighbours than the 
        # rest, this way of storing the neighbour list mat be inneficient
        self.neighlist = np.full((self.nspins, nneigh_max), -1, dtype="intc")

        # Fill the list
        for j_spin in range(self.nspins):
            neighs = self.network.neighbours_out(j_spin)
            self.neighlist[j_spin][0:neighs.size] = neighs

        return 

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

            Returns
            -------
                naccept : int
                    Number of accepted proposals.

        """
        # Raise an error if the given temperature is not valid
        if T <= 0:
            raise ValueError("The temperature must be greater than zero.")

        # Calculate the value of the thermodynamic beta
        # for the given temperature
        beta = 1./T

        # Call the time evolution C function
        self.spins, naccept = cevolve.evolve_nofieldGlauber(
                self.spins, self.neighlist, self.nneigh, beta, nsteps)

        return naccept


    # TODO: improve docs
    def slow_thermalization(self, T_final, T_ini=4., T_step=0.1,
                            steps_per_T=1000):
        """Thermalize the system slowly to the given temperature.

        Parameters
       ----------
            T_final : float
                Final temperature.

            T_ini : float
                Initial temperature in the thermalization process.

            T_step : float
                Difference between consecutive temperatures.

            steps_per_T : int
                Number of time steps per temperature.

            T_step : float
                Difference between consecutive temperature steps.

        """
        # Create the temperature vector
        Ts = np.arange(T_ini, T_final, T_step)

        # Temperature loop
        for T in Ts:
            self.evolve(steps_per_T, T)

        return

    # Measures
    # ========================================
    def magnetization(self):
        """Calculate the magnetization.

        """
        return np.mean(self.spins)

    # TODO: check how slow is this function an use cython if required
    def hamiltonian(self):
        """Calculate the value of the hamiltonian.

        All couplings are 1 and no external field is considered.

        """
#        hamilt = 0.
#        for pair in self.neighpairs:
#            hamilt += -self.spins[pair[0]]*self.spins[pair[1]]
        # Call the C function
        hamilt = cevolve.hamiltonian(self.spins, self.neighpairs)
        
        return hamilt


    # Graphical representation
    # ========================================

    def plot(self, figsize=3):
        """Plot the system configuration. 

        This method only works if the lattice is bidimensional.

        """

        fig, ax = plt.subplots(figsize=(figsize, figsize))
        im = ax.imshow(self.latt, cmap="Greys", vmin=-1, vmax=+1,
                interpolation=None)
        return fig


    def animate(self, T, nframes=1000, steps_per_frame=1, frame_interval=300):
        """Animate the evolution of the lattice. 

        This method only works if the lattice is bidimensional.

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


    # I/O
    # ========================================

    def saveconf(self, fname=None):
        """Save spin configuration of the system to a text file.

        Parameters
        ----------
            fname : string
                Name of the output file. Defaults to spinconfN@.dat
                with @ the number of spins in the lattice.

        """
        if fname == None:
            fname = "spinconfN{0}.dat".format(self.nspins)

        np.savetxt(filename, latt.spins, fmt="%d")
        return

    def loadconf(self, fname):
        """Load spin configuration of the system from a text file.

        Parameters
        ----------
            fname : string
                Name of the input file. 

        """
        inconf = np.loadtxt(filename).astype("intc")

        # If the number of spins does not match the number of
        # spins in the lattice, return an error
        if inconf.size != self.nspins:
            raise ValueError(
                    "The number of spins in the input file does not match the "
                    "number of spins in the lattice.")

        latt.spins = inconf

        return

class IsingCoupling(Ising):
    """Ising model class with arbitrary couplings.

    The coupling are taken from the adjacency matrix of the network
    attribute. The adjacency matrix array must be of type "float64".

    """
    def evolve(self, nsteps, T):
        """Evolve the lattice in nsteps.

        Accept arbitrary couplings.

           Parameters
           ----------
               nsteps : int
                   Number of Monte Carlo steps.

               T : float
                   Temperature of the system in units of the inverse
                   of the Boltzmann's constant times the coupling
                   constant times the particle's spin squared.

            Returns
            -------
                naccept : int
                    Number of accepted proposals.

        """
        # Raise an error if the given temperature is not valid
        if T <= 0:
            raise ValueError("The temperature must be greater than zero.")

        # Calculate the value of the thermodynamic beta
        # for the given temperature
        beta = 1./T

        # Call the time evolution C function
        self.spins, naccept = cevolve.evolve_nofieldCouplGlauber(
                self.spins, self.neighlist, self.network.adjmatrix, self.nneigh, 
                beta, nsteps)

        return naccept

class Regular(Ising):
    """Regular Ising model lattice with periodic boundary conditions.

    """
    def __init__(self, shape, seed=None):
        """ Init method.

        Parameters
        ----------
            shape : int list
                Lattice shape.

            seed : int
                Random number generator seed.

        """
        self.nspins = np.prod(shape)

        # Create the regular network
        network = networks.Lattice(
                shape, pbc=True, weighted=False, directed=False)

        Ising.__init__(self, self.nspins, network=network, shape=shape, seed=seed)


class Ising2D(Regular):
    """2D Ising model lattice with periodic boundary conditions.

    """
    def __init__(self, nrows, ncols, seed=None):
        # Store parameters
        shape = tuple((nrows, ncols))

        Regular.__init__(self, shape, seed=seed)
        

#    @latt.setter
#    def latt(self, value):
#        np.reshape(self.spins, self.shape()) = value


#    def neighbours(self):
#        """Create neighbour list for the lattice.
#
#        Each spin has its the closest ones (up, down, left, right)
#        has neighbours, with periodic boundary condiions taked 
#        into account.
#
#        Returns
#        -------
#            neighlist : int array
#                2D array with the neighbours of each cell.
#                neighbourlist[i,j] is the index in self.latt
#                of the j-th neighbour of the i-th cell (with
#                i also its index in self.latt). The j values
#                correspond, in increasing order, to the
#                top, right, bottom and left neighbours.
#
#            nneigh: int 
#                Number of neighbours of each cell.
#        
#        """
#        # Number of neighbours per spin.
#        nneigh = 4 
#        neighlist = np.zeros((self.nspins, nneigh), dtype="intc")
#
#        for j_spin in range(self.nspins):
#            # Calculate the row and column of the spin.
#            [row, col] = np.unravel_index(j_spin, self.shape())
#
#            # Store the neighbours' indexes.
#
#            # Use the wrap parameter to implement the periodic boundaries.
#            neighlist[j_spin][0] = np.ravel_multi_index(
#                    [row-1, col], self.shape(), mode="wrap") # Top neighbour
#            neighlist[j_spin][1] = np.ravel_multi_index(
#                    [row, col+1], self.shape(), mode="wrap") # Right neighbour
#            neighlist[j_spin][2] = np.ravel_multi_index(
#                    [row+1, col], self.shape(), mode="wrap") # Bottom neighbour
#            neighlist[j_spin][3] = np.ravel_multi_index(
#                    [row, col-1], self.shape(), mode="wrap") # Left neighbour
#            
#        return neighlist, nneigh


    # Exact results
    # ========================================

    # Critical temperature
    Tcrit = 2.26918531421

    # Correlation length critical exponent
    corrlen_exp = 1.
    # Magnetization critical exponent
    mag_exp = 1./8
    # Magnetic critical exponent
    magsuscept_exp = 7./4
    # Specific susceptibility exponent
    specificheat_exp = 0.
    
    @classmethod
    def mag_exact(cls, T):
        """Return the expected magnetization using the exact solution.

        Parameters
        ----------
            T : float (or float array)
                Temperature.

        Returns
        -------
            mag : float
                Magnetization.

        """
        # Convert T to a numpy array
        # TODO: there must be a cleaner way to do this
        T_arr = T*np.ones(1)

        # Create array to store the solution
        mag = np.zeros(T_arr.size, dtype=float)

        # Find the indexes of the temperatures below the critical one
        undercrit_idxs = np.argwhere(T_arr < cls.Tcrit)
        # Calculate the magnetization for the undercritical temperatures
        beta = 1./T_arr[undercrit_idxs]
        aux = np.power(np.sinh(2.*beta), -4)
        mag[undercrit_idxs] = np.power(1. - aux, 1./8.)

        # If there is only one temperature, convert the array to a scalar
        if T_arr.size == 1:
            mag = np.asscalar(mag)
        return mag
