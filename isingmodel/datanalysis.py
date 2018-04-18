import numpy as np
import random 
from matplotlib import pyplot as plt
from matplotlib import animation

from memoize import memoize


class Results(object):
    # TODO: improve docs
    def __init__(self, shape=None, fname=None, nsigma=1.):
        """Blalbalba

        Parameters
        ----------
            shape : int 2-tuple
                Shape of the lattice whose measures are stored.

            fname : string
                Name of a text file to be imported.

            nsigma : float
                The error in a measured magnitudes will be nsigma
                times the standard deviation.

        """
        # Store parameters
        self.nsigma = nsigma

        if shape != None:
            self.shape = list(shape)
        else:
            self.shape = None

        # If the filename is provided, read the data from there
        if fname != None:
            self.readtxt(fname)

        else:
            # Store parameters
            if self.shape == None:
                raise ValueError("Lattice shape not given.")

            # Initialize results lists
            self.Ts = list()
            self.mags = list()
            self.mag2s = list()
            self.mag4s = list()
            self.corrmags = list()
            self.hamilts = list()
            self.hamilt2s = list()
            self.hamilt4s = list()
            self.corrhamilts = list()
            self.nmeasures = list()
            self.acceptprobs = list()
            self.measureintervals = list()

        # Calculate the numer of spins
        self.nspins = np.prod(self.shape)

    # TODO: complete docs
    # TODO: check if T has been already measured and average
    # with the previous data in that case
    def measure(self, T, nmeasures, latt, measureinterval=1):
        """Measure blablbalba

        """
        # Check if lattice shape is the expected one
        if self.shape != latt.shape():
            raise ValueError(
                    "The lattice shape does not match the Results object one.")
            
        # Store parameters
        self.Ts.append(T)
        self.nmeasures.append(nmeasures)
        self.measureintervals.append(measureinterval)

        # Initialize variables
        mag_last = 0. # Magnetization in the last measure
        hamilt_last = 0. # Hamiltonian in the last measure
        mag_sum = 0.
        mag2_sum = 0.
        mag4_sum = 0.
        corrmag_sum = 0.
        hamilt_sum = 0.
        hamilt2_sum = 0.
        hamilt4_sum = 0.
        corrhamilt_sum = 0.
        naccept = 0

        # Start measure loop
        for measure_idx in range(nmeasures):
            # Evolve
            naccept += latt.evolve(measureinterval, T) 

            # Measure
            mag = latt.magnetization()
            mag2 = mag*mag
            hamilt = latt.hamiltonian()
            hamilt2 = hamilt*hamilt

            mag_sum += np.abs(mag)
            mag2_sum += mag2
            mag4_sum += mag2*mag2
            corrmag_sum += mag*mag_last
            hamilt_sum += hamilt
            hamilt2_sum += hamilt2
            hamilt4_sum += hamilt2*hamilt2
            corrhamilt_sum += hamilt*hamilt_last

            # Store last measure
            mag_last = mag
            hamilt_last = hamilt

        # Store measures and calculate means
        self.mags.append(mag_sum/nmeasures)
        self.mag2s.append(mag2_sum/nmeasures)
        self.mag4s.append(mag4_sum/nmeasures)
        self.corrmags.append(corrmag_sum/(nmeasures - 1))
        self.hamilts.append(hamilt_sum/nmeasures)
        self.hamilt2s.append(hamilt2_sum/nmeasures)
        self.hamilt4s.append(hamilt4_sum/nmeasures)
        self.corrhamilts.append(corrhamilt_sum/(nmeasures - 1))
        self.acceptprobs.append(
                float(naccept)/(nmeasures*measureinterval*latt.nspins))

        return


    @memoize
    def L(self):
        """Return characteristic size of the system.

        """
        return np.power(np.prod(self.shape), 1./self.shape.size)


    # I/O 
    # ==============================
    # TODO: add the data instead of overwriting it and check if the shape
    # of the imported file is the same as the object attribute
    def readtxt(self, filename):
        """Read data from file.

        """
        filedata = np.loadtxt(filename).T

        self.Ts = filedata[0].tolist()
        self.mags = filedata[1].tolist()
        self.mag2s = filedata[2].tolist()
        self.mag4s = filedata[3].tolist()
        self.corrmags = filedata[4].tolist()
        self.hamilts = filedata[5].tolist()
        self.hamilt2s = filedata[6].tolist()
        self.hamilt4s = filedata[7].tolist()
        self.corrhamilts = filedata[8].tolist()
        self.acceptprobs = filedata[9].tolist()
        self.nmeasures = filedata[10].tolist()
        self.measureintervals = filedata[11].tolist()

        # Read additional parameters from footer
        with open(filename, "r") as f:
            lines = f.readlines()
            self.shape = map(int, lines[-1].split()[2].split(","))

        return

    def savetxt(self, fname=None):
        """Save data to file.

        Parameters
        ----------
            fname : string
                Name of the output file. Its default value is
                "isingR{0}C{1}.dat" with {0} the number of rows 
                in the lattice and {1} the number of columns.

        """
        if fname == None:
            fname = "isingR{0}C{1}.dat".format(self.shape[0], self.shape[1])

        headerstring = (
                "Temperature\t "
                "Mean mag.\t Mag. 2nd moment\t Mag. 4nd moment\t "
                "Mag. time corr.\t " 
                "Mean hamilt.\t Hamilt. 2nd moment\t Hamilt. 4nd moment\t "
                "Hamilt. time corr.\t " 
                "Acceptance probability\t N measures\t Measure interval")
        footerstring = "Shape: {0},{1}".format(self.shape[0], self.shape[1])

        np.savetxt(
                fname,
                np.vstack((
                        self.Ts, self.mags, self.mag2s, self.mag4s,
                        self.corrmags, self.hamilts, self.hamilt2s, 
                        self.hamilt4s, self.corrhamilts, self.acceptprobs,
                        self.nmeasures, self.measureintervals)).T,
                header=headerstring, footer=footerstring)
        return

        
    # Physical magnitudes
    # ========================================
    def mag_err(self):
        """Calculate the magnetization error.

        """
        # Calculate correlation time
        corrtime = corr_time(
                self.mags, self.mag2s, self.corrmags, self.nmeasures)
        return self.nsigma*samplemean_error(
                self.mags, self.mag2s, corrtime, self.nmeasures)

    def mag2_err(self):
        """Calculate the error of the squared magnetization mean.

        """
        # Calculate correlation time. We are making the assumtion
        # that the correlation time of mag2 is the same as mag.
        corrtime = corr_time(
                self.mags, self.mag2s, self.corrmags, self.nmeasures)
        return self.nsigma*samplemean_error(
                self.mag2s, self.mag4s, corrtime, self.nmeasures)

    def hamilt_err(self):
        """Calculate the Hamiltonian error.

        """
        # Calculate correlation time
        corrtime = corr_time(
                self.hamilts, self.hamilt2s, self.corrhamilts, self.nmeasures)
        return self.nsigma*samplemean_error(
                self.hamilts, self.hamilt2s, corrtime, self.nmeasures)

    def hamilt2_err(self):
        """Calculate the error of the squared Hamiltonian mean.

        """
        # Calculate correlation time. We are making the assumtion
        # that the correlation time of hamilt2 is the same as hamilt's.
        corrtime = corr_time(
                self.hamilts, self.hamilt2s, self.corrhamilts, self.nmeasures)
        return self.nsigma*samplemean_error(
                self.hamilt2s, self.hamilt4s, corrtime, self.nmeasures)

    def magsuscept(self):
        """Calculate the magnetic susceptibility.

        """
        # Store data to numpy arrays
        Ts_arr = np.array(self.Ts)
        return self.nspins/Ts_arr*samplevariance(
                self.mags, self.mag2s, self.nmeasures)

    def magsuscept_err(self):
        """Calculate the magnetic susceptibility error.

        """
        # Store data to numpy arrays
        Ts_arr = np.array(self.Ts)

        return self.nspins/Ts_arr*np.sqrt(
                np.power(self.mag2_err(), 2)
                + 4.*np.power(self.mags*self.mag_err(), 2))

    def specificheat(self):
        """Calculate the specific heat per spin of the lattice.

        """
        # Store data to numpy arrays
        Ts_arr = np.array(self.Ts)

        return 1./(self.nspins*np.power(Ts_arr, 2))*samplevariance(
                self.hamilts, self.hamilt2s, self.nmeasures)

    def specificheat_err(self):
        """Calculate the specific heat per spin error.

        """
        # Store data to numpy arrays
        Ts_arr = np.array(self.Ts)

        return 1./(self.nspins*np.power(Ts_arr, 2))*np.sqrt(
                np.power(self.hamilt2_err(), 2)
                + 4.*np.power(self.hamilts*self.hamilt_err(), 2))

                
    def binderratio(self):
        """Calculate the Binder ratio or fourth order cumulant.

        """
        return (1. - self.mag4s/(3.*np.power(self.mag2s, 2)))

    

# Statistical functions
# ===================================
def variance(mean, momnt2):
    """Calculate the sample variance.

    Parameters
    ----------
        mean : float (scalar or array)
            Mean value. 

        momnt2 : float (scalar or array)
            Second raw moment (mean of the square).

    Returns
    -------
        variance : float (scalar or array)

    """
    momnt2_arr = np.array(momnt2)
    return momnt2_arr - np.power(mean, 2)

def samplevariance(mean, momnt2, nmeasure):
    """Calculate the sample variance.

    Parameters
    ----------
        mean : float (scalar or array)
            Mean value. 

        momnt2 : float (scalar or array)
            Second raw moment (mean of the square).

    Returns
    -------
        variance : float (scalar or array)

    """
    nmeasure_arr = np.array(nmeasure)
    return nmeasure_arr/(nmeasure_arr - 1.)*variance(mean, momnt2)


# TODO: improve docs
# TODO: ensure the units are right
def corr_time(mean, momnt2, corr, nmeasures):
    """Estimate the correlation time in a Markov chain (with rejection).

    Estimates the correlation time using the mean value
    of the product in consecutive steps and the variance
    (it is assumed that the autocorrelation decays
    exponentially).
    
    Parameters
    ----------
        mean : float (scalar or array)
            Mean of the magnitued.

        momnt2 : float (scalar or array)
            Second moment of the magnitude.

        corr : float (scalar or array)
            Mean value of the product of the magnitude in
            consecutive measures.

        nmeasures: int (scalar or array)
            Number of measures.

    Returns
    -------
        corr_time : float (scalar or array)
            Estimated correlation time.
        
    """
    # Calculate the variance
    var = samplevariance(mean, momnt2, nmeasures)

    # Ensure the data is stored in arrays
    var_arr = var*np.ones(1)
    corr_arr = corr*np.ones(1)
    mean_arr = mean*np.ones(1)

    # Find the indexes where the variance is not zero
    nonzero_idxs = np.argwhere(var_arr != 0)

    # Initialize to -1
    corr_norm = np.full(corr_arr.shape, -1., dtype=float)

    # Calculate the normalized autocorrelation
    corr_norm[nonzero_idxs] = (
            (corr_arr[nonzero_idxs] - np.power(mean_arr[nonzero_idxs], 2))
            /var_arr[nonzero_idxs])
    return corr_norm/(1. - corr_norm)


def samplemean_error(mean, momnt2, corrtime, nmeasures):
    """Calculate the sample mean error in rejection with repetition.

    Parameters
    ----------
        mean : float (scalar or array)
            Sample mean of the calculated magnitued.

        momnt2 : float (scalar or array)
            Sample second raw moment of the magnitude.

        corrtime : float (scalar or array)
            Correlation time of the magnitude.

        nmeasures: int (scalar or array)
            Number of measures.

    Returns
    -------
        error : float (scalar or array)
        
    """
    # Calculate the variance
    var = samplevariance(mean, momnt2, nmeasures)

    # If the variance is zero, the error is directly zero.
    # If we use the formula in those cases a zero division is
    # done, so we have to treat the zero values separately.

    # Ensure the data is stored in arrays
    mean_arr = mean*np.ones(1)
    var_arr = var*np.ones(1)
    corrtime_arr = corrtime*np.ones(1)
    nmeasures_arr = nmeasures*np.ones(1)
    
    # Create array for the results
    error = np.zeros(var_arr.size, dtype=float)

    # Find the array indexes with nonzero variance and calculate
    # the error in those cases
    nonzero_idxs = np.argwhere(var_arr != 0)
    error[nonzero_idxs] = np.sqrt(
            var_arr[nonzero_idxs]/nmeasures_arr[nonzero_idxs]*(
            2.*corrtime_arr[nonzero_idxs] + 1.))

    # If the array size is one, convert it to a scalar
    if error.size == 1:
        error = np.asscalar(error)

    return error


# Auxiliar functions
# ========================================

# TODO: make the function check and treat properly measures with 
# the same T but different measure intervals
# TODO: there must be a way of doing this cleanly
def mergeresults(results_list):
    """Merge several results objects into one.

    Be careful: right now the function does not treat properly 
    measures with the same temperature but different 
    measure intervals.

    Parameters
    ----------
        results_list : isingmodel.Results list
            List with the Results objects to be merged. All the 
            Results objects in the list bust have the same shape.
            in the list must have the 

    Returns
    -------
        merged: :py:class:`isingmodel.Results` object
            Results object with all the data from results list.

    """
    # Check that all the list elements have the same shape
    shape = results_list[0].shape
    for results in results_list:
        if results.shape != shape:
            raise ValueError(
                    "All the elements in the list must have the same shape")

    # Create the object where all the data will be stored
    merged = Results(shape=shape)

    # Loop over the results objects
    for results in results_list:
        # Loop over the measured temperatures
        for T_idx, T in enumerate(results.Ts):
            # If T is already in merged, average the results
            if T in merged.Ts:
                idx = merged.Ts.index(T) 

                merged.nmeasures[idx] += results.nmeasures[T_idx]
                merged.mags[idx] = np.average(
                        [merged.mags[idx], results.mags[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
                merged.mag2s[idx] = np.average(
                        [merged.mag2s[idx], results.mag2s[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
                merged.mag4s[idx] = np.average(
                        [merged.mag4s[idx], results.mag4s[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
                merged.corrmags[idx] = np.average(
                        [merged.corrmags[idx], results.corrmags[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
                merged.hamilts[idx] = np.average(
                        [merged.hamilts[idx], results.hamilts[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
                merged.hamilt2s[idx] = np.average(
                        [merged.hamilt2s[idx], results.hamilt2s[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
                merged.hamilt4s[idx] = np.average(
                        [merged.hamilt4s[idx], results.hamilt4s[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
                merged.corrhamilts[idx] = np.average(
                        [merged.corrhamilts[idx], results.corrhamilts[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
                merged.acceptprobs[idx] = np.average(
                        [merged.acceptprobs[idx], results.acceptprobs[T_idx]],
                        weights=[merged.nmeasures[idx], 
                                results.nmeasures[T_idx]])
            # Else, create a new entry for the measures
            else:
                merged.Ts.append(
                        results.Ts[T_idx])
                merged.nmeasures.append(
                        results.nmeasures[T_idx])
                merged.measureintervals.append(
                        results.measureintervals[T_idx])
                merged.mags.append(
                        results.mags[T_idx])
                merged.mag2s.append(
                        results.mag2s[T_idx])
                merged.mag4s.append(
                        results.mag4s[T_idx])
                merged.corrmags.append(
                        results.corrmags[T_idx])
                merged.hamilts.append(
                        results.hamilts[T_idx])
                merged.hamilt2s.append(
                        results.hamilt2s[T_idx])
                merged.hamilt4s.append(
                        results.hamilt4s[T_idx])
                merged.corrhamilts.append(
                        results.corrhamilts[T_idx])
                merged.acceptprobs.append(
                        results.acceptprobs[T_idx])
    return merged
