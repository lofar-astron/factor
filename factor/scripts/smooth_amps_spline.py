#! /usr/bin/env python
"""
Script to smooth and normalize amplitude solutions
Reinout van Weeren, April 2016
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import numpy
import os
import lofar.parmdb
import math
import shutil
import multiprocessing
import matplotlib.pyplot as plt
from scipy.interpolate import LSQUnivariateSpline
import sys
import scipy.ndimage
import astropy.convolution
import matplotlib as mpl


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
    scatter = numpy.median(abs(shifted_vec - datavector))
    return scatter


def findscatter_time(dataarray):
    scattervec = []
    for freq in range(0,len(dataarray[:,0])):
      #print 'findscatter_time', freq
      scatter = findscatter(dataarray[freq,:])
      scattervec.append(scatter)
    return numpy.median(scattervec)


def findscatter_freq(dataarray):
    scattervec = []
    for time in range(0,len(dataarray[0,:])):
      #print 'findscatter_freq', time
      scatter = findscatter(dataarray[:,time])
      scattervec.append(scatter)
    return numpy.median(scattervec)


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

    # work in log-sapce
    amp_orig_ext = numpy.copy(amp)
    amp = numpy.log10(amp)
    weights = (0.*numpy.copy(amp)) + 1 # initialize weights to 1

    # filter bad data and determine average scatter of amplitudes
    idx = numpy.where(amp != 0.0) # log10(1.0) = 0.0
    #print idx, 'idx'
    if numpy.any(idx): # so we do not have an empty array
        scatter = findscatter(amp[idx])
        # remove some really bad stuff, by putting weights to zero.
        idxbadi1 = numpy.where(amp > (numpy.median(amp) + (35.*std(amp))))
        weights[idxbadi1] = 1e-10 # small value, zero generates NaN in spline
        idxbadi2 = numpy.where(amp < (numpy.median(amp) - (35.*std(amp))))
        weights[idxbadi2] = 1e-10  # small value, zero generates NaN in spline
    else:
        scatter = 0.02 # just that we have a value to prevent crashes in case all amplitudes are 1.0
        #print 'No valid data for found for this anntenna: ', antenna


    # make the noisevec
    if numpy.any(idx): # at least 1 good data point
        if len(amp[idx]) > 30:  # so at least 30/3 = 10 good data points
            # create noise vector
            noisevec = findnoisevec(amp)
        else:
            noisevec = (numpy.copy(amp) * 0.) + 1.0 # just make constant noise, if we have too little datapoints
    else:
        noisevec = (numpy.copy(amp) * 0.) + 1.0 # just make constant noise, if we have too little datapoints


    #print scatter, antenna

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

    #print 'knots CL', knotvec
    spl2 = LSQUnivariateSpline(timevec, amp, knotvec, w=weights, k=splineorder)

    # now find bad data devatiating from the fit 30 x scatter
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
    #idxnodata = numpy.where(amp_orig_ext == 1.0)
    #amp[idxnodata] = 0.0 # to avoid problem with amplitudes that are 1.0
    amp = 10**amp

    amp_clean = amp[ndata:ndata + ndata]

    idxbad = numpy.where(amp_clean != amp_orig)
    n_knots = numpy.int(numpy.ceil(numpy.float(len(knotvec))/3.)) # approxmiate, just for plot

    # return cleaned amplitudes, model, scatter, number of knots, indices of replaced outliers
    return amp_clean, 10**(model[ndata:ndata + ndata]), noisevec[ndata:ndata + ndata], scatter, n_knots, idxbad, weights[ndata:ndata + ndata]


def pad_2Darray(a, width, mode):
    pad_shape = (a.shape[0]*3, a.shape[1]*3)
    pad_a = np.zeros(pad_shape)

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

    # take the log
    amp = numpy.log10(amp)

    # create median filtered array
    amp_median = scipy.ndimage.median_filter(amp, (3,5)) # so a bit more smoothing along the time-axis

    # find scatter
    idxgood = numpy.where(amp != 0.0)
    if numpy.any(idxgood):
        scatter_freq = findscatter_freq(amp)
        scatter_time = findscatter_time(amp)
    else:
        scatter_freq  = 0.02 # just asign some value
        scatter_time =  0.02 # just asign some value
        #print 'scatter (freq,time)', scatter_freq, scatter_time, antenna, pol

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


def main(instrument_name, instrument_name_smoothed, normalize=True, plotting=False):
    if type(normalize) is str:
        if normalize.lower() == 'true':
            normalize = True
        else:
            normalize = False
    if type(plotting) is str:
        if plotting.lower() == 'true':
            plotting = True
        else:
            plotting = False

    gain = 'Gain'

    pdb = lofar.parmdb.parmdb(instrument_name)
    parms = pdb.getValuesGrid('*')

    key_names = parms.keys()
    nchans = len(parms[key_names[0]]['freqs'])

    # determine the number of polarizations in parmdb (2 or 4)
    if any(gain+':0:1:' in s for s in key_names):
        pol_list = ['0:0', '1:1', '0:1', '1:0']
    else:
        pol_list = ['0:0', '1:1']

    times = numpy.copy(sorted( parms[key_names[0]]['times']))
    freqs = numpy.copy(sorted( parms[key_names[0]]['freqs']))/1e6 # get this in MHz

    # times not used at the moment, I assume the time axis for a parmdb is regular and does not contain gaps
    times = (times - numpy.min(times))/24. #so we get an axis in hrs

    # Get station names
    antenna_list = set([s.split(':')[-1] for s in pdb.getNames()])

    # for plotting
    Nr = int(numpy.ceil(numpy.sqrt(len(antenna_list))))
    Nc = int(numpy.ceil(numpy.float(len(antenna_list))/Nr))
    if plotting:
        mpl.rc('font',size =6 )
        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95 )
        fa, axa = plt.subplots(Nr, Nc, sharex=True, sharey=True, figsize=(16,12))
        axsa = axa.reshape((Nr*Nc,1))
        if nchans > 5: # 2D filter
            Nr = len(antenna_list)
            Nc = int(4)
            fa2, axa2 = plt.subplots(Nr, Nc, sharex=True, sharey=True, figsize=(8,108),)
            axsa2 = axa2.reshape((Nr*Nc,1))

    for pol in pol_list:
        for istat,antenna in enumerate(sorted(antenna_list)[::-1]):
            channel_parms_real = [parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan]
                for chan in range(nchans)]
            channel_parms_imag = [parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan]
                for chan in range(nchans)]

            # some plotting setup
            if len(channel_parms_real[0]) > 500:
                fmt = ','
            else:
                fmt = 'o'
            ls='none'

            channel_amp_orig = [numpy.sqrt(channel_parms_real[chan]**2 +
                channel_parms_imag[chan]**2) for chan in range(nchans)]

            # now find the bad data
            pool = multiprocessing.Pool()
            results = pool.map(spline1D, channel_amp_orig)
            pool.close()
            pool.join()

            for chan, (amp_cleaned, model, noisevec, scatter, n_knots, idxbad, weights) in enumerate(results):
                # put back the results
                phase = numpy.arctan2(channel_parms_imag[chan], channel_parms_real[chan]**2)
                parms[gain + ':' + pol + ':Real:' + antenna]['values'][:, chan] = numpy.copy(amp_cleaned*numpy.cos(phase))
                parms[gain + ':' + pol + ':Imag:' + antenna]['values'][:, chan] = numpy.copy(amp_cleaned*numpy.sin(phase))

                if pol in pol_list[0]:
                    cc = 'blue'
                    ccf = 'orange'
                else:
                    cc = 'green'
                    ccf= 'red'

                timevec = numpy.arange(0,len(channel_amp_orig[chan]))

                # only plot one channel, just to verify code works
                if plotting and chan == nchans-1: # plot last channel
                    axsa[istat][0].plot(timevec, amp_cleaned, marker=fmt, ls=ls,
                        markersize=0.1*len(amp_cleaned), c=cc,mec=cc)
                    axsa[istat][0].plot(timevec,noisevec, c=cc, lw=0.75, ls='--')

                    if pol in pol_list[0]:
                        axsa[istat][0].annotate('scatter=' +'{:.2g}'.format(scatter),
                            xy=(0.5,0.15), color=cc,textcoords='axes fraction')
                        axsa[istat][0].annotate('#knots=' +'{:d}'.format(n_knots),
                            xy=(0.01,0.15), color=cc,textcoords='axes fraction') # we divded by three beucase we mirrored the array
                    else:
                        axsa[istat][0].annotate('scatter=' +'{:.2g}'.format(scatter),
                            xy=(0.5,0.02), color=cc, textcoords='axes fraction')
                        axsa[istat][0].annotate('#knots=' +'{:d}'.format(n_knots),
                            xy=(0.01,0.02), color=cc,textcoords='axes fraction')

                    if numpy.any(idxbad):
                        axsa[istat][0].plot(timevec[idxbad],channel_amp_orig[chan][idxbad],
                            marker='o', c=ccf, ls=ls, markersize=4)

                    idxbadi = numpy.where(weights < 1.0)
                    if numpy.any(idxbadi):
                        axsa[istat][0].plot(timevec[idxbadi],channel_amp_orig[chan][idxbadi],
                            marker='o', c='black', ls=ls, markersize=4, mec='black')

                    axsa[istat][0].plot(timevec, model, c=ccf, lw=1.0)
                    axsa[istat][0].set_title(antenna)
                    axsa[istat][0].set_ylim(-0.3, 2)
                    axsa[istat][0].set_xlim(0, max(timevec))

            if nchans > 5: # Do 2D smooth
                channel_parms_real = [parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan]
                    for chan in range(nchans)]
                channel_parms_imag = [parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan]
                    for chan in range(nchans)]
                channel_parms_real =  numpy.asarray(channel_parms_real)
                channel_parms_imag =  numpy.asarray(channel_parms_imag)
                channel_amp_orig = [numpy.sqrt(channel_parms_real[chan]**2 +
                    channel_parms_imag[chan]**2) for chan in range(nchans)]
                amp_orig = numpy.sqrt(channel_parms_real[:]**2 + channel_parms_imag[:]**2)
                phase    = numpy.arctan2(channel_parms_imag[:], channel_parms_real[:]**2)

                amp_cleaned, amp_median, baddata = median2Dampfilter(numpy.copy(amp_orig))

                for chan in range(nchans):
                    # put back the results
                    parms[gain + ':' + pol + ':Real:' + antenna]['values'][:, chan] = numpy.copy((amp_cleaned[chan,:])*numpy.cos(phase[chan,:]))
                    parms[gain + ':' + pol + ':Imag:' + antenna]['values'][:, chan] = numpy.copy((amp_cleaned[chan,:])*numpy.sin(phase[chan,:]))

                if plotting:
                    axsa2[4*istat][0].imshow(numpy.transpose(amp_orig),
                        interpolation='none',origin='lower',clim=(0.5, 1.5),aspect='auto')
                    axsa2[4*istat][0].set_xlabel('freq')
                    axsa2[4*istat][0].set_ylabel('time')
                    axsa2[4*istat][0].set_title('Original' + '    ' + antenna)

                    axsa2[4*istat+1][0].imshow(numpy.transpose(amp_median),
                        interpolation='none',origin='lower',aspect='auto', clim=(0.5,1.5))
                    axsa2[4*istat+1][0].set_xlabel('freq')
                    axsa2[4*istat+1][0].set_ylabel('time')
                    axsa2[4*istat+1][0].set_title('2D median model')

                    axsa2[4*istat+2][0].imshow(numpy.transpose(numpy.abs(amp_orig-amp_median)),
                        interpolation='none',origin='lower',clim=(0.0, 0.3),aspect='auto')
                    axsa2[4*istat+2][0].set_xlabel('freq')
                    axsa2[4*istat+2][0].set_ylabel('time')
                    axsa2[4*istat+2][0].set_title('abs(Residual)')

                    axsa2[4*istat+3][0].imshow(numpy.transpose(baddata),
                        interpolation='none',origin='lower',clim=(0.0, 2.0),
                        aspect='auto', cmap='gnuplot')
                    axsa2[4*istat+3][0].set_xlabel('freq')
                    axsa2[4*istat+3][0].set_ylabel('time')
                    axsa2[4*istat+3][0].set_title('Replaced solutions')


    if plotting:
        fa.savefig('1Dsmooth.png', dpi=100)
        if nchans > 5: # make 2D plot
            fa2.tight_layout()
            fa2.savefig('2Dsmooth.png')
        plt.show()

    # Normalize the amplitude solutions to a mean of one across all channels
    if normalize:
        # First find the normalization factor
        amplist = []
        for chan in range(nchans):
            for pol in ['0:0','1:1']:  # hard code here in case the data contains 0:1 and 1:0
                for antenna in antenna_list:
                    real = numpy.copy(parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan])
                    imag = numpy.copy(parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan])
                    amp  = numpy.copy(numpy.sqrt(real**2 + imag**2))
                    amplist.append(amp)
        norm_factor = 1.0/(numpy.mean(amplist))
        print "smooth_amps_spline.py: Normalization-Factor is:", norm_factor

        # Now do the normalization
        for chan in range(nchans):
            for pol in pol_list:
                for antenna in antenna_list:
                    real = numpy.copy(parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan])
                    imag = numpy.copy(parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan])
                    phase = numpy.arctan2(imag, real)
                    amp  = numpy.copy(numpy.sqrt(real**2 + imag**2))

                    # Clip extremely low amplitude solutions to prevent very high
                    # amplitudes in the corrected data
                    low_ind = numpy.where(amp < 0.2)
                    amp[low_ind] = 0.2

                    parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan] = numpy.copy(amp *
                        numpy.cos(phase) * norm_factor)
                    parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan] = numpy.copy(amp *
                        numpy.sin(phase) * norm_factor)

    if os.path.exists(instrument_name_smoothed):
        shutil.rmtree(instrument_name_smoothed)
    pdbnew = lofar.parmdb.parmdb(instrument_name_smoothed, create=True)
    pdbnew.addValues(parms)
    pdbnew.flush()


if __name__ == '__main__':
    descriptiontext = "Smooth and normalize amplitude solutions.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('instrument_name', help='name of the instrument parmdb to smooth')
    parser.add_argument('instrument_name_smoothed', help='name of the output parmdb')
    parser.add_argument('normalize', help='normalize?')
    parser.add_argument('plotting', help='make plots?')
    args = parser.parse_args()

    main(args.instrument_name, args.instrument_name_smoothed, args.normalize, args.plotting)
