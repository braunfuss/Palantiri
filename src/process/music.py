#   Adapted by Andreas Steinberg 2019 for teleseismic BP
#   based on 2013 Russell Haley doamusic tool

from __future__ import absolute_import, division, print_function
import scipy as sp
import numpy as np
import itertools
from scipy import linalg
from scipy import misc
from scipy import pi



def _spectrum(
    metric,
    antennas,
    out,
    thlo,thstep,thsz,
    phlo,phstep,phsz
):
    # Lower-level spectrum calculator with preprocessed arguments and
    # pass-by-reference output array, for easier implementation with
    # cython and being farmed out to multiple processes. (The problem is
    # embarassingly parallel.
    assert out.shape == (thsz,phsz)
    for i in range(thsz):
        th = thlo + i*thstep
        for j in range(phsz):
            ph = phlo + j*phstep
            out[i,j] = _pmusic(metric,antennas,th,ph)


def sph2cart(sph):
    """
    Convert one or more spherical coordinates to cartesian.

    Parameters
    ----------
    sph : shape 3 or Nx3 numpy.ndarray
        Spherical coordinates of the form (r,theta,phi), where theta is the
        inclination angle from the Z axis and phi is the azimuth angle from the
        X axis.

    Returns
    -------
    cart : shape 3 or Nx3 numpy.ndarray
        Cartesian coordinates of the form (x,y,z)
    """
    sph = np.atleast_2d(sph).T
    x = sph[0] * np.sin(sph[1]) * np.cos(sph[2])
    y = sph[0] * np.sin(sph[1]) * np.sin(sph[2])
    z = sph[0] * sp.cos(sph[1])
    cart = np.squeeze(np.array((x,y,z)).T)
    return cart

def cart2sph(cart):
    """
    Convert one or more cartesian coordinates to spherical.

    Parameters
    ----------
    cart : shape 3 or Nx3 numpy.ndarray
        Cartesian coordinates of the form (x,y,z)

    Returns
    -------
    sph : shape 3 or Nx3 numpy.ndarray
        Spherical coordinates of the form (r,theta,phi), where theta is the
        inclination angle from the Z axis and phi is the azimuth angle from the
        X axis.
    """
    cart = np.atleast_2d(cart).T
    r = np.sqrt(np.sum(cart**2,axis=0))
    th = np.arctan2( np.sqrt(np.sum(cart[:2]**2,axis=0)), cart[2] )
    ph = np.arctan2(cart[1],cart[0])
    cart = np.squeeze(np.array((r,th,ph)).T)
    return cart


def eigsort(eigresult):
    """
    Sort the output of scipy.linalg.eig() in terms of
    eignevalue magnitude
    """
    ix = sp.argsort(abs(eigresult[0]))
    return ( eigresult[0][ix], eigresult[1][:,ix] )


class Estimator:
    """
    A class to carry state for estimating direction of arrival (DOA)
    with the Multiple SIgnal Classification algorithm.
    """
    def __init__(
        self,
        antennas,
        covariance,
        field_of_view=((0,pi),(-pi,pi)),
        nsignals=None
    ):
        """
        Set up an Estimator, for making pseudospectrum plots and finding
        directions of arrival.

        Parameters
        ----------
        antennas : sequence of of N sequences of 3, or Nx3 numpy array
            Describes the relative physical positions of the array elements.
            Units are wavelength.

        covariance : NxN numpy array
           The joint second moment matrix.  The ij entry is the expectation of
           the product of the conjugated signal from the ith element with the
           signal from the jth element. Equal to mean[sample.H * sample].

        field_of_view : ((th_lo,th_hi),(ph_lo,ph_hi))
            Restrict the domain of DoA searches and pseudospectrum renderings.
            This is helpful to avoid aliasing from 2-dimensional arrays, or to
            avoid projection distortion in the pseudospectrum image by
            rendering only a subsection of the sphere.  Theta and phi are the
            inclination angle from the Z axis and azimuth angle from the X
            axis.

        nsignals : integer < N
            The number of incident signals, if known. Otherwise, we will try to
            estimate this from the magnitudes of the eigenvalues of the
            covariance matrix.
        """
        # Accept and validate antennas.
        self.antennas = np.array(antennas).astype(complex) * (2*pi)
        self.numel = self.antennas.shape[0]
        assert self.antennas.shape[1] == 3      # we are operating in R3
        # Accept and validate covariance.
        self.covar = np.array(covariance)
        assert self.covar.shape == (self.numel,self.numel)
        # Unpack field of view
        self.thlo,self.thhi = field_of_view[0][0],field_of_view[0][1]
        self.phlo,self.phhi = field_of_view[1][0],field_of_view[1][1]

        #Get the sorted eigenstructure
        self.eigval, self.eigvec = eigsort(linalg.eig(covariance))

        # Try to guess the number of incident signals, if unspecified
        if nsignals:
            assert nsignals < self.numel
            self.nsignals = nsignals
            self.noisedim = self.numel - nsignals
        else:
            shaped = abs(self.eigval)
            self.noisedim = sp.diff(shaped).argmax() + 1
            self.nsignals = self.numel - self.noisedim

        #slice the noise space
        self.noisespace = self.eigvec[:,:self.noisedim]
        self.sigspace = self.eigvec[:,self.noisedim:]
        #DEBUG print("Noise space dimension: {}".format(self.noisespace.shape))

        # Calculate the noise metric used to evaluate pmusic, to avoid
        # repetition
        self.metric = sp.atleast_2d(
                    self.noisespace.dot( self.noisespace.T.conj() )
                 ).astype(complex)

    def eigplot():
        """
        Plot the eigenvalues on a logarithmic chart, for human appraisal
        of the number of incident signals (and any other interesting
        properties of the data.
        """
        pass

    def spectrum(self,shape,method=_spectrum):
        """
        Generate a MUSIC pseudospectrum on the estimator's domain. The result
        is a theta_sz x phi_sz real numpy.ndarray. The domain is a closed
        interval, like linspace.

        Parameters
        ----------
        shape : (theta_sz, phi_sz)
            Specify the dimensions of the result

        method : callable
            Choose between the python or cython low-level implementations.
            Used to check correctness.
        """
        # Wraps either _spectrum or _music.spectrum and provides parallel
        # evaluation.

        # Extract shape argument
        theta_sz,phi_sz = shape

        # precalculate static arguments as comlpex double and prepare output
        # array
        result = np.empty((theta_sz,phi_sz))

        # step sizes
        thstep = (self.thhi-self.thlo)/(theta_sz-1)
        phstep = (self.phhi-self.phlo)/(phi_sz-1)

        method(
           self.metric,
           self.antennas,
           result,
           self.thlo,thstep,theta_sz,
           self.phlo,phstep,phi_sz
        )
        return result

    def doasearch(self,max_iterations=2**8,tol=pi/2**16):
        """
        Find directions of arrival within specified tolerance.
        """
        out = []
        iterations = 0
        # Try to find twice as many signals, so we get the aliases too.
        while iterations < max_iterations:
            # track the number of iterations
            iterations += 1
            # Pick a random starting point on the sphere by taking the theta
            # and phi coordinates of a random sample from the normal
            # distribution (which is spherically symmetric).
            thstart,phstart = cart2sph(sp.randn(3))[1:]
            phstart = 0 + sp.rand()*pi
            maxval,maxpoint = _music.hillclimb(
                self.metric,
                self.antennas,
                thstart,phstart,
                2**(-80)
            )
            for val,point in out:
                # Check distance from ones we've already found.
                diff = sp.array(point) - sp.array(maxpoint)
                if sp.sqrt(sp.dot(diff,diff)) < 2*tol:
                    break
            else:
                # It's not one of the ones we've already found.
                out.append((maxval,maxpoint))
        # Sort the signals by descending pmuisc amplitude
        out.sort(key=lambda x: x[0], reverse=True)
        # Filter the points outside the domain.
        out = [ point for val,point in out
            if (
                point[0] < self.thhi and
                point[0] > self.thlo and
                point[1] < self.phhi and
                point[1] > self.phlo
            )
        ]
        # Take the signals biggest (which should hopefully discard aliases).
        out = out[:self.nsignals]
        #DEBUG: print("Iterations: {}".format(iterations))
        return out




def aoa2prop(th,ph):
    """Propogation vectors from angles of arrival. Result is Nx3."""
    return np.array((-np.sin(th)*np.cos(ph), #x
                     -np.sin(th)*np.sin(ph), #y
                     -np.cos(th)             #z
                   )).T # we want Nx3 not 3xN

def aoa2prop_scalar(th,ph):
    """Slightly faster for single arguments, 'cause numpy scalars are slow."""
    return np.array((-math.sin(th)*math.cos(ph), #x
                     -math.sin(th)*math.sin(ph), #y
                     -math.cos(th)               #z
                   ))


def covar(samples):
    """
    Calculate the covariance matrix as used by Estimator.

    This is not the same as the Octave/Matlab function cov(), but is instead
    equal to Mean [ sample.H * sample ], where sample is a single sample.
    I.E., it is actually the second moment matrix.

    Parameters
    ----------
    samples : K x Numel or Numel x 0 complex ndarray
        Complex samples for each of Numel antennas sampled at K times.

    Returns
    -------
        return : Numel x Numel complex ndarray
            Second moment matrix for complex random vector samples.  Used by
            Estimator.
    """
    samples = sp.asmatrix(samples)
    return ( (samples.H * samples) / samples.shape[0] )

def _pmusic(metric,antennas,theta,phi):
    steer = sp.exp( 1j*antennas.dot(-aoa2prop_scalar(theta,phi)) )
    return 1.0 / steer.conj().dot(metric).dot(steer).real
