#! /usr/bin/env python
"""
Script to smooth and normalize amplitude solutions
Reinout van Weeren, April 2016
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import numpy
import os
from losoto.h5parm import h5parm
import math
import shutil
import multiprocessing
from scipy.interpolate import LSQUnivariateSpline, interp1d, interp2d
import sys
import scipy.ndimage
import astropy.convolution


def std(inputData, Zero=False, axis=None, dtype=None):
    """
    Robust estimator of the standard deviation of a data set.

    Based on the robust_sigma function from the AstroIDL User's Library.

    .. versionchanged:: 1.0.3
        Added the 'axis' and 'dtype' keywords to make this function more
        compatible with numpy.std()
    """
    epsilon = 1.0e-20
    if axis is not None:
        fnc = lambda x: std(x, dtype=dtype)
        sigma = numpy.apply_along_axis(fnc, axis, inputData)
    else:
        data = inputData.ravel()
        if type(data).__name__ == "MaskedArray":
            data = data.compressed()
        if dtype is not None:
            data = data.astype(dtype)

        if Zero:
            data0 = 0.0
        else:
            data0 = numpy.median(data)
        maxAbsDev = numpy.median(numpy.abs(data-data0)) / 0.6745
        if maxAbsDev < epsilon:
            maxAbsDev = (numpy.abs(data-data0)).mean() / 0.8000
        if maxAbsDev < epsilon:
            sigma = 0.0
            return sigma

        u = (data-data0) / 6.0 / maxAbsDev
        u2 = u**2.0
        good = numpy.where( u2 <= 1.0 )
        good = good[0]
        if len(good) < 3:
            print "WARNING:  Distribution is too strange to compute standard deviation"
            sigma = -1.0
            return sigma

        numerator = ((data[good]-data0)**2.0 * (1.0-u2[good])**2.0).sum()
        nElements = (data.ravel()).shape[0]
        denominator = ((1.0-u2[good])*(1.0-5.0*u2[good])).sum()
        sigma = nElements*numerator / (denominator*(denominator-1.0))
        if sigma > 0:
            sigma = math.sqrt(sigma)
        else:
            sigma = 0.0

    return sigma

def findscatter(datavector):
    shifted_vec = numpy.roll(datavector, 1)
    #scatter = sum(abs(shifted_vec - datavector))/numpy.float(len(datavector))
    scatter = numpy.nanmedian(abs(shifted_vec - datavector))
    return scatter


def findscatter_time(dataarray):
    scattervec = []
    for freq in range(0,len(dataarray[:,0])):
      #print 'findscatter_time', freq
      scatter = findscatter(dataarray[freq,:])
      scattervec.append(scatter)
    return numpy.nanmedian(scattervec)


def findscatter_freq(dataarray):
    scattervec = []
    for time in range(0,len(dataarray[0,:])):
      #print 'findscatter_freq', time
      scatter = findscatter(dataarray[:,time])
      scattervec.append(scatter)
    return numpy.nanmedian(scattervec)


def findnoisevec(datavector):
    shifted_vec = numpy.roll(datavector, 1)
    scatter_vec = (abs(shifted_vec - datavector))
    #scatter_vec = medfilt(scatter_vec,21)
    scatter_vec = scipy.ndimage.filters.median_filter(scatter_vec,9, mode='mirror')

    # now smooth
    gauss = astropy.convolution.Gaussian1DKernel(stddev=4.0)
    scatter_vec = astropy.convolution.convolve(scatter_vec,gauss , boundary='extend')

    # normalize scatter_vec
    scatter_vec = scatter_vec/numpy.mean(scatter_vec)

    return scatter_vec


def spline1D(amp_orig):
    # to compute knot points
    f = lambda m, n: [i*n//m + n//(2*m) for i in range(m)]

    if amp_orig is None:
        return None, None, None, None, None, None, None

    # expand array and mirror full array around edges
    ndata = len(amp_orig)
    amp = numpy.zeros(ndata+2*ndata)
    amp[ndata:ndata+ndata] = amp_orig

    for i in range(0, ndata):
        # Mirror at left edge.
        idx = min(ndata-1, ndata-i)
        amp[i] = amp_orig[idx]
        # Mirror at right edge
        idx = max(0, ndata-2-i)
        amp[ndata+ndata+i] = amp_orig[idx]

    # Find flagged values
    flagged = numpy.where(amp == 1.0)

    # work in log-sapce
    amp_orig_ext = numpy.copy(amp)
    amp = numpy.log10(amp)
    weights = (0.*numpy.copy(amp)) + 1 # initialize weights to 1

    # filter bad data and determine average scatter of amplitudes
    scatter = findscatter(amp)
    # remove some really bad stuff, by putting weights to zero.
    idxbadi1 = numpy.where(amp > (numpy.median(amp) + (35.*std(amp))))
    weights[idxbadi1] = 1e-10 # small value, zero generates NaN in spline
    idxbadi2 = numpy.where(amp < (numpy.median(amp) - (35.*std(amp))))
    weights[idxbadi2] = 1e-10  # small value, zero generates NaN in spline

    # Set weights for flagged values
    weights[flagged] = 1e-10  # small value, zero generates NaN in spline

    # make the noisevec
    if len(amp) > 30:  # so at least 30/3 = 10 good data points
        # create noise vector
        noisevec = findnoisevec(amp)
    else:
        noisevec = (numpy.copy(amp) * 0.) + 1.0 # just make constant noise, if we have too little datapoints

    if scatter < 0.005:
        #Interior knots t must satisfy Schoenberg-Whitney conditions
        scatter = 0.005 # otherwise we fit more parameters than we have data points
    knotfactor = 0.5e3*scatter  # normalize based on trial and error

    #print scatter, antenna, len(amp), knotfactor

    timevec = numpy.arange(0,len(amp))
    knotvec = f(numpy.int(len(amp)/knotfactor),len(amp))
    #print antenna, 'knots', knotvec, noisevec[knotvec]

    #print 'knots OR', knotvec

    # simple optimization knot selection for vectors that have at least 30 data points
    # based on the noisevector
    # removes even numbered knots if the noise is high
    knotvec_copy = numpy.copy(knotvec) # otherwise tcopy is updated as well
    if len(timevec) > 30 and len(knotvec) > 2:
        for counter, knot in enumerate(knotvec_copy):
            #print counter, knot, noisevec[knot]
            if (counter % 2 == 0) and noisevec[knot] > 1.5: # even index and large noise
                knotvec.remove(knot)
              #print 'Removing knot because of local increase in noise'

    #print antenna, 'cleaned knots', knotvec, noisevec[knotvec]

    # asign midpoint if not enough data points/20
    if len (knotvec) < 3: # because we are working with a 3x larger mirrored array
        knotvec = [numpy.int(len(timevec)*0.25),numpy.int(len(timevec)/2),numpy.int(len(timevec)*0.75)]
        #print 'extending to', knotvec

    splineorder =  5 #  default
    if len(knotvec) == 3 and scatter > 0.1:
        splineorder = 3 # reduce order, data is  bad
        if scatter > 0.2:
            splineorder = 1 # very bad data
    spl2 = LSQUnivariateSpline(timevec, amp, knotvec, w=weights, k=splineorder)

    # now find bad data devatiating from the fit 15 x scatter
    residual = numpy.abs(spl2(timevec)-amp)
    idx      = numpy.where(residual > 15.*scatter)

    # second iteration
    if numpy.any(idx):
        ampcopy = numpy.copy(amp)
        ampcopy[idx] = spl2(timevec[idx]) # replace bad amplitudes by model
        spl2 = LSQUnivariateSpline(timevec, ampcopy, knotvec,  w=weights, k=splineorder)

    residual = numpy.abs(spl2(timevec)-amp)
    idx      = numpy.where(residual > 8.*scatter)

    # third iteration
    if numpy.any(idx):
        ampcopy = numpy.copy(amp)
        ampcopy[idx] = spl2(timevec[idx]) # replace bad amplitudes by model
        spl2 = LSQUnivariateSpline(timevec, ampcopy, knotvec,  w=weights, k=splineorder)

    # again look at residual, go back to original amp again, find deviating data > 3x scatter
    residual = numpy.abs(spl2(timevec)-amp)
    idx      = numpy.where(residual > 3.*scatter)
    # replace the bad data with model
    model    =spl2(timevec)

    #if len(idx) != 0:
    amp[idx] = model[idx]

    # go out of log-space
    idxnodata = numpy.where(numpy.logical_or(amp > 1.0, amp < -10.0))
    amp[idxnodata] = 0.0
    amp[flagged] = 0.0
    model[flagged] = 0.0
    amp = 10**amp

    amp_clean = amp[ndata:ndata + ndata]

    idxbad = numpy.where(amp_clean != amp_orig)
    n_knots = numpy.int(numpy.ceil(numpy.float(len(knotvec))/3.)) # approxmiate, just for plot

    # return cleaned amplitudes, model, scatter, number of knots, indices of replaced outliers
    return amp_clean, 10**(model[ndata:ndata + ndata]), noisevec[ndata:ndata + ndata], scatter, n_knots, idxbad, weights[ndata:ndata + ndata]


def pad_2Darray(a, width, mode):
    pad_shape = (a.shape[0]*3, a.shape[1]*3)
    pad_a = numpy.zeros(pad_shape)

    # center
    pad_a[a.shape[0]:2*a.shape[0], a.shape[1]:2*a.shape[1]] = a

    # four corners
    pad_a[0:a.shape[0], 0:a.shape[1]] = a[::-1, ::-1]
    pad_a[0:a.shape[0], 2*a.shape[1]:3*a.shape[1]] = a[::-1, ::-1]
    pad_a[2*a.shape[0]:3*a.shape[0], 2*a.shape[1]:3*a.shape[1]] = a[::-1, ::-1]
    pad_a[2*a.shape[0]:3*a.shape[0], 0:a.shape[1]] = a[::-1, ::-1]

    # middle edges
    pad_a[0:a.shape[0], a.shape[1]:2*a.shape[1]] = a[:, ::-1]
    pad_a[a.shape[0]:2*a.shape[0], 2*a.shape[1]:3*a.shape[1]] = a[::-1, :]
    pad_a[2*a.shape[0]:3*a.shape[0], a.shape[1]:2*a.shape[1]] = a[:, ::-1]
    pad_a[a.shape[0]:2*a.shape[0], 0:a.shape[1]] = a[::-1, :]

    return pad_a


def median2Dampfilter(amp_orig):
    try:
        from numpy import pad
    except ImportError:
        pad = pad_2Darray

    orinal_size = numpy.shape(amp_orig)
    # padd array by reflection around axis

    amp = pad(amp_orig, ((numpy.shape(amp_orig)[0],numpy.shape(amp_orig)[0]),
        (numpy.shape(amp_orig)[1],numpy.shape(amp_orig)[1])), mode='reflect')
    flagged = numpy.where(numpy.logical_or(amp == 1.0, amp == 0.0))

    # take the log
    amp = numpy.log10(amp)

    # Set flagged values to NaN
    amp[flagged] = numpy.nan

    # create median filtered array, ignoring NaNs
    amp_median = scipy.ndimage.filters.generic_filter(amp, numpy.nanmedian, (3,5)) # so a bit more smoothing along the time-axis

    # find scatter
    scatter_freq = findscatter_freq(amp)
    scatter_time = findscatter_time(amp)
    # print 'scatter (freq,time)', scatter_freq, scatter_time

    scatter = 0.5*(scatter_freq+scatter_time) # average x-y scatter

    # find bad data
    idxbad = numpy.where((numpy.abs(amp - amp_median)) > scatter*3.)
    baddata = numpy.copy(amp)*0.0
    baddata[idxbad] = 1.0

    # replace the bad data points
    amp_cleaned = numpy.copy(amp)
    amp_cleaned[idxbad] = amp_median[idxbad]

    # raise to the power
    amp = 10**amp
    amp_median = 10**amp_median
    amp_cleaned = 10**amp_cleaned

    #back to original size
    amp_median = amp_median[orinal_size[0]:2*orinal_size[0],orinal_size[1]:2*orinal_size[1]]
    baddata   = baddata[orinal_size[0]:2*orinal_size[0],orinal_size[1]:2*orinal_size[1]]
    amp_cleaned = amp_cleaned[orinal_size[0]:2*orinal_size[0],orinal_size[1]:2*orinal_size[1]]

    return amp_cleaned, amp_median, baddata


def main(instrument_name, normalize=True):
    if type(normalize) is str:
        if normalize.lower() == 'true':
            normalize = True
        else:
            normalize = False

    h = h5parm(instrument_name)
    solset = H.getSolset('sol000')
    soltab = solset.getSoltab('amplitude000')
    parms = soltab.val[:]

    # Check for NaNs and zeros. If found, set to 1
    ntimes = len(soltab.time[:])
    nfreqs = len(soltab.freq[:])
    initial_flagged_indx = numpy.where(numpy.logical_or(numpy.isnan(parms), parms == 0.0))
    initial_unflagged_indx = numpy.where(numpy.logical_and(~numpy.isnan(parms), parms != 0.0))
    parms[initial_flagged_indx] = 1.0

    # Get station names
    antenna_list = soltab.ant[:]
    axis_names = soltab.getAxesNames()
    for pol in len(soltab.pol[:]):
        for istat, antenna in enumerate(sorted(antenna_list)[::-1]):
            soltab.setSelection(ant=[antenna])
            channel_amp_orig = soltab.val[:] # ['time', 'freq', 'ant', 'dir', 'pol']

            # Find flagged solutions and set to 1.0
            channel_amp_interp = []
            for chan in range(nchans):
                unflagged_times = numpy.where(channel_amp_orig[:, chan, 0, 0, pol] != 1.0)
                flagged_times = numpy.where(channel_amp_orig[:, chan, 0, 0, pol] == 1.0)
                if numpy.any(flagged_times):
                    channel_amp_orig[flagged_times, chan, 0, 0, pol] = 1.0
                if numpy.any(unflagged_times):
                    channel_amp_interp.append(channel_amp_orig[:, chan, 0, 0, pol])
                else:
                    channel_amp_interp.append(None)

            # now find the bad data
            if ntimes > 5:
                pool = multiprocessing.Pool()
                results = pool.map(spline1D, channel_amp_interp)
                pool.close()
                pool.join()

                for chan, (amp_cleaned, model, noisevec, scatter, n_knots, idxbad, weights) in enumerate(results):
                    # put back the results
                    if amp_cleaned is None:
                        amp_cleaned = channel_amp_orig[:, chan, 0, 0, pol]
                    parms[:, chan, istat, 0, pol] = numpy.copy(amp_cleaned)

            if nchans > 5: # Do 2D smooth
                channel_amp_orig = parms[:, :, istat, 0, pol]

                # Smooth
                channel_amp_interp = []
                unflagged_sols = numpy.where(channel_amp_orig != 1.0)
                if numpy.any(unflagged_sols):
                    # Set flagged solutions to 1.0
                    flagged_sols = numpy.where(numpy.logical_or(channel_amp_orig == 1.0, channel_amp_orig <= 0.0))
                    channel_amp_orig[flagged_sols] = 1.0

                    # Filter
                    amp_cleaned, amp_median, baddata = median2Dampfilter(channel_amp_orig.transpose([1, 0]))
                    amp_cleaned = amp_cleaned.transpose([1, 0])
                    amp_cleaned[flagged_sols] = 1.0
                    parms[:, :, istat, 0, pol] = numpy.copy(amp_cleaned)

    # Normalize the amplitude solutions to a mean of one across all channels
    if normalize:
        # First find the normalization factor from unflagged solutions
        amplist = []
        norm_factor = 1.0/(numpy.nanmean(parms[initial_unflagged_indx]))
        print "smooth_amps_spline.py: Normalization-Factor is:", norm_factor
        parms *= norm_factor

        # Clip extremely low amplitude solutions to prevent very high
        # amplitudes in the corrected data
        unflagged = numpy.where(~numpy.isnan(parms))
        low_ind = numpy.where(parms[unflagged] < 0.2)
        amp[unflagged][parms] = 0.2

    # Make sure flagged solutions are still flagged
    parms[initial_flagged_indx] = numpy.nan

    # Write the results to the output soltab
    st = solset.makeSoltab('amplitude', 'amplitude001', axesNames=axis_names,
                           axesVals=soltab.getAxisValues(), vals=parms,
                           weights=np.ones(parms.shape))

if __name__ == '__main__':
    descriptiontext = "Smooth and normalize amplitude solutions.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('instrument_name', help='name of the instrument h5parm to smooth')
    parser.add_argument('normalize', help='normalize?')
    args = parser.parse_args()

    main(args.instrument_name, args.normalize)
