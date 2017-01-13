#! /usr/bin/env python

# Stat definitions

import numpy as np
import numpy.ma as ma

def chi2(obs,err):
  obs = obs/np.median(obs)
  core = ((obs - 1.) / err)**2.
  return core.sum()/len(obs-1)

def chi2_master_ref(obs,err,comp,comp_err):
  #obs = obs/np.median(obs)
  #comp = comp/np.median(comp)
  core = ((obs - comp) / np.sqrt(err**2+comp_err**2))**2.
  return core.sum()/len(obs-1)

def robust(obs,err):
  core = abs((obs-np.median(obs))/err)
  return core.sum()/len(obs-1)

def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default
    """
    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m

def coloumn_median(B,A):
  maxlen = max(len(x) for x in B)
  if len(A) == maxlen:
    ans = A
  else:
    C = np.array([l+[np.nan]*(maxlen-len(l)) for l in B])
    dat = np.ma.masked_array(C,np.isnan(C))
    ans = np.median(dat,axis=0)
  return ans

def coloumn_weighted_mean(B,err,N_ref_stars):
  data = {}
  weights_N = {}
  results = []
  max_len = 0
  
  weights = [[] for _ in range(N_ref_stars)]
  for i in range(len(err)):
    err[i] = np.array(err[i])
    variance = err[i]*err[i]
    weights[i] = 1./variance
    weights[i] = weights[i]/weights[i].sum()

  for alist in B:
	  length = len(alist)
	  max_len = length if (length > max_len) else max_len
	  for i in range(length):
		  data.setdefault(i, []).append(alist[i])

  for alist in weights:
	  length = len(alist)
	  max_len = length if (length > max_len) else max_len
	  for i in range(length):
		  weights_N.setdefault(i, []).append(alist[i])
		  
  for i in range(max_len):
	  vals = np.array(data[i])
	  vals_w = np.array(weights_N[i])
	  results.append(np.average(vals, weights=vals_w))
  
  return np.array(results)

def sigma_clip(data, sig=3, iters=1, cenfunc=np.median, varfunc=np.var,
               maout=False):
    """ Perform sigma-clipping on the provided data.

    This performs the sigma clipping algorithm - i.e. the data will be iterated
    over, each time rejecting points that are more than a specified number of
    standard deviations discrepant.

    .. note::
        `scipy.stats.sigmaclip` provides a subset of the functionality in this
        function.

    Parameters
    ----------
    data : array-like
        The data to be sigma-clipped (any shape).
    sig : float
        The number of standard deviations (*not* variances) to use as the
        clipping limit.
    iters : int or None
        The number of iterations to perform clipping for, or None to clip until
        convergence is achieved (i.e. continue until the last iteration clips
        nothing).
    cenfunc : callable
        The technique to compute the center for the clipping. Must be a
        callable that takes in a 1D data array and outputs the central value.
        Defaults to the median.
    varfunc : callable
        The technique to compute the variance about the center. Must be a
        callable that takes in a 1D data array and outputs the width estimator
        that will be interpreted as a variance. Defaults to the variance.
    maout : bool or 'copy'
        If True, a masked array will be returned. If the special string
        'inplace', the masked array will contain the same array as `data`,
        otherwise the array data will be copied.

    Returns
    -------
    filtereddata : `numpy.ndarray` or `numpy.masked.MaskedArray`
        If `maout` is True, this is a masked array with a shape matching the
        input that is masked where the algorithm has rejected those values.
        Otherwise, a 1D array of values including only those that are not
        clipped.
    mask : boolean array
        Only present if `maout` is False. A boolean array with a shape matching
        the input `data` that is False for rejected values and True for all
        others.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    only the points that are within 2 *sample* standard deviation from the
    median::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> data,mask = sigma_clip(randvar, 2, 1)

    This will clipping on a similar distribution, but for 3 sigma relative to
    the sample *mean*, will clip until converged, and produces a
    `numpy.masked.MaskedArray`::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> maskedarr = sigma_clip(randvar, 3, None, mean, maout=True)

    """

    data = np.array(data, copy=False)
    oldshape = data.shape
    data = data.ravel()

    mask = np.ones(data.size, bool)
    if iters is None:
        i = -1
        lastrej = sum(mask) + 1
        while(sum(mask) != lastrej):
            i += 1
            lastrej = sum(mask)
            do = data - cenfunc(data[mask])
            mask = do * do <= varfunc(data[mask]) * sig ** 2
        iters = i + 1
        #TODO: ?print iters to the log if iters was None?
    else:
        for i in range(iters):
            do = data - cenfunc(data[mask])
            mask = do * do <= varfunc(data[mask]) * sig ** 2

    if maout:
        return np.ma.MaskedArray(data, ~mask, copy=maout != 'inplace')
    else:
        return data[mask], mask.reshape(oldshape)
