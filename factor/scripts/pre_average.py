#! /usr/bin/env python
"""
Script to pre-average data using a sliding Gaussian kernel on the weights
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import glob
import sys
import os
import itertools
import pickle
from scipy.ndimage.filters import gaussian_filter1d as gfilter
import pyrap.tables as pt
import lofar.parmdb
from astropy.stats import median_absolute_deviation


def main(ms_file, parmdb_file, input_colname, output_colname, output_weights_colname,
    pre_average=True, minutes_per_block=10.0, baseline_file=None, verbose=True):
    """
    Pre-average data using a sliding Gaussian kernel on the weights
    """
    if not pre_average:
        # Just copy input columns to output
        ms = pt.table(ms_file, readonly=False, ack=False)
        if output_colname not in ms.colnames():
            desc = ms.getcoldesc(input_colname)
            desc['name'] = output_colname
            ms.addcols(desc)
        if output_weights_colname not in ms.colnames():
            desc = ms.getcoldesc('WEIGHT_SPECTRUM')
            desc['name'] = output_weights_colname
            ms.addcols(desc)
        data = ms.getcol(input_colname)
        weights = ms.getcol('WEIGHT_SPECTRUM')
        ms.putcol(output_colname, data)
        ms.putcol(output_weights_colname, weights)
        ms.flush()
        ms.close()
    else:
        # Do the BL averaging
        if baseline_file is None:
            if verbose:
                print('Calculating baseline lengths...')
            baseline_dict = get_baseline_lengths(ms_file)
        elif os.path.exists(baseline_file):
            f = open('baseline_file', 'r')
            baseline_dict = pickle.load(f)
            f.close()
        else:
            print('Cannot find baseline_file. Exiting...')
            sys.exit(1)

        # Iterate through time chunks and find the lowest ionfactor
        tab = pt.table(ms_file, ack=False)
        start_time = tab[0]['TIME']
        end_time = tab[-1]['TIME']
        remaining_time = end_time - start_time # seconds
        tab.close()
        ionfactors = []
        t_delta = minutes_per_block * 60.0 # seconds
        t1 = 0.0
        if verbose:
            print('Determining ionfactors...')
        while remaining_time > 0.0:
            if remaining_time < 1.5 * t_delta:
                # If remaining time is too short, just include it all in this chunk
                t_delta = remaining_time + 10.0
            remaining_time -= t_delta

            # Find ionfactor for this period
            ionfactors.append(find_ionfactor(parmdb_file, baseline_dict, t1+start_time,
                t1+start_time+t_delta))
            if verbose:
                print('    ionfactor (for timerange {0}-{1} sec) = {2}'.format(t1,
                    t1+t_delta, ionfactors[-1]))
            t1 += t_delta

        # Do pre-averaging using lowest ionfactor
        ionfactor_min = min(ionfactors)
        if verbose:
            print('Using ionfactor = {}'.format(ionfactor_min))
            print('Averaging...')
        BLavg(ms_file, baseline_dict, input_colname, output_colname,
            output_weights_colname, ionfactor_min)


def get_baseline_lengths(ms_file):
    """
    Returns dict of baseline lengths in km for all baselines in input dataset
    """
    t = pt.table(ms_file, ack=False)
    anttab = pt.table(ms_file+'::ANTENNA', ack=False)
    antnames = anttab.getcol('NAME')
    anttab.close()
    ant1 = t.getcol('ANTENNA1')
    ant2 = t.getcol('ANTENNA2')
    all_uvw = t.getcol('UVW')
    baseline_dict = {}
    for ant in itertools.product(set(ant1), set(ant2)):
        if ant[0] >= ant[1]:
            continue
        sel1 = np.where(ant1 == ant[0])[0]
        sel2 = np.where(ant2 == ant[1])[0]
        sel = sorted(list(frozenset(sel1).intersection(sel2)))
        uvw = all_uvw[sel, :]
        uvw_dist = np.sqrt(uvw[:, 0]**2 + uvw[:, 1]**2 + uvw[:, 2]**2)
        baseline_dict['{0}'.format(ant[0])] = antnames[ant[0]]
        baseline_dict['{0}'.format(ant[1])] = antnames[ant[1]]
        baseline_dict['{0}-{1}'.format(ant[0], ant[1])] = np.mean(uvw_dist) / 1.e3
    t.close()

    return baseline_dict


def find_ionfactor(parmdb_file, baseline_dict, t1, t2):
    """
    Finds ionospheric scaling factor
    """
    pdb_in = lofar.parmdb.parmdb(parmdb_file)
    parms = pdb_in.getValuesGrid('*')

    # Filter any stations not in both the instrument table and the ms
    stations_pbd = set([s.split(':')[-1] for s in pdb_in.getNames()])
    stations_ms = set([s for s in baseline_dict.itervalues() if type(s) is str])
    stations = sorted(list(stations_pbd.intersection(stations_ms)))

    # Select long baselines only, as they will set the ionfactor scaling
    ant1 = []
    ant2 = []
    dist = []
    min_length = 10.0
    for k, v in baseline_dict.iteritems():
        if type(v) is not str and '-' in k:
            if v > min_length:
                s1 = k.split('-')[0]
                s2 = k.split('-')[1]
                s1_name = baseline_dict[s1]
                s2_name = baseline_dict[s2]
                if s1_name in stations and s2_name in stations:
                    ant1.append(s1_name)
                    ant2.append(s2_name)
                    dist.append(v)

    # Find correlation times
    target_rms_rad = 0.2
    rmstimes = []
    freq = None
    for a1, a2, d in zip(ant1, ant2, dist):
        if freq is None:
            freq = np.copy(parms['Gain:0:0:Phase:{}'.format(a1)]['freqs'])[0]
            times = np.copy(parms['Gain:0:0:Phase:{}'.format(a1)]['times'])
            time_ind = np.where((times >= t1) & (times < t2))[0]
            timepersolution = np.copy(parms['Gain:0:0:Phase:{}'.format(a1)]['timewidths'])[0]
        ph1 = np.copy(parms['Gain:0:0:Phase:{}'.format(a1)]['values'])[time_ind]
        ph2 = np.copy(parms['Gain:0:0:Phase:{}'.format(a2)]['values'])[time_ind]

        rmstime = None
        ph = unwrap_fft(ph2 - ph1)

        step = 1
        for i in range(1, len(ph)/2, step):
            p1 = ph[i:]
            p2 = ph[:-i]
            rms = np.linalg.norm(p1-p2) / np.sqrt(len(p1))
            mad = median_absolute_deviation(p1-p2)
            mean = np.mean(p1-p2)
            if rms + mean > target_rms_rad:
                rmstime = i
                break
        if rmstime is None:
            rmstime = len(ph)/2
        rmstimes.append(rmstime)

    # Find the mean ionfactor assuming that the correlation time goes as
    # t_corr ~ 1/sqrt(BL). The ionfactor is defined in BLavg() as:
    #
    #     ionfactor = (t_corr / 30.0 sec) / ( np.sqrt((25.0 / dist_km)) * (freq_hz / 60.e6) )
    #
    ionfactor = np.mean(np.array(rmstimes) / 30.0 / (np.sqrt(25.0 / np.array(dist))
        * freq / 60.0e6)) * timepersolution

    return ionfactor


def BLavg(msfile, baseline_dict, input_colname, output_colname, output_weights_colname,
    ionfactor, clobber=True):
    """
    Averages data using a sliding Gaussian kernel on the weights
    """
    if not os.path.exists(msfile):
        print("Cannot find MS file.")
        sys.exit(1)

    # open input/output MS
    ms = pt.table(msfile, readonly=False, ack=False)
    freqtab = pt.table(msfile + '::SPECTRAL_WINDOW', ack=False)
    freq = freqtab.getcol('REF_FREQUENCY')
    freqtab.close()
    wav = 299792458. / freq
    timepersample = ms.getcell('INTERVAL',0)
    all_time = ms.getcol('TIME_CENTROID')

    ant1 = ms.getcol('ANTENNA1')
    ant2 = ms.getcol('ANTENNA2')
    all_data = ms.getcol(input_colname)
    all_weights = ms.getcol('WEIGHT_SPECTRUM')
    all_flags = ms.getcol('FLAG')

    # iteration on baseline combination
    for ant in itertools.product(set(ant1), set(ant2)):

        if ant[0] >= ant[1]:
            continue
        sel1 = np.where(ant1 == ant[0])[0]
        sel2 = np.where(ant2 == ant[1])[0]
        sel = sorted(list(frozenset(sel1).intersection(sel2)))

        # compute the FWHM
        dist = baseline_dict['{0}-{1}'.format(ant[0], ant[1])]
        stddev = 30.0 * ionfactor * np.sqrt((25.0 / dist)) * (freq / 60.e6) # in sec
        stddev = stddev/timepersample # in samples

        #    Multiply every element of the data by the weights, convolve both
        #    the scaled data and the weights, and then divide the convolved data
        #    by the convolved weights (translating flagged data into weight=0).
        #    That's basically the equivalent of a running weighted average with
        #    a Gaussian window function.

        # get weights
        flags = all_flags[sel,:,:]
        weights = all_weights[sel,:,:]*~flags # set flagged data weight to 0
        # get data
        data = all_data[sel,:,:]*weights

        # smear weighted data and weights
        dataR = gfilter(np.real(data), stddev, axis=0)#, truncate=4.)
        dataI = gfilter(np.imag(data), stddev, axis=0)#, truncate=4.)
        weights = gfilter(weights, stddev, axis=0)#, truncate=4.)

        # re-create data
        # NOTE: both data and/or weight might be 0
        d = (dataR + 1j * dataI)
        weights[np.where(data == 0)] = np.nan # prevent 0/0 exception
        d = d/weights # when weights == 0 or nan -> data is nan
        # if data is nan, put data to 0 (anyway they must be flagged)
        all_data[sel,:,:] = np.nan_to_num(d)
        all_weights[sel,:,:] = np.nan_to_num(weights)

    # Add the output columns if needed
    if output_colname not in ms.colnames():
        desc = ms.getcoldesc(input_colname)
        desc['name'] = output_colname
        ms.addcols(desc)
    if output_weights_colname not in ms.colnames():
        desc = ms.getcoldesc('WEIGHT_SPECTRUM')
        desc['name'] = output_weights_colname
        ms.addcols(desc)

    ms.putcol(output_colname, all_data)
    ms.putcol(output_weights_colname, all_weights)
    ms.close()


def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    import numpy as np
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))

    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]


def unwrap_fft(phase, iterations=3):
    """
    Unwrap phase using Fourier techniques.

    For details, see:
    Marvin A. Schofield & Yimei Zhu, Optics Letters, 28, 14 (2003)

    Keyword arguments:
    phase -- array of phase solutions
    iterations -- number of iterations to perform
    """
    puRadius=lambda x : np.roll( np.roll(
          np.add.outer( np.arange(-x.shape[0]/2+1,x.shape[0]/2+1)**2.0,
                        np.arange(-x.shape[1]/2+1,x.shape[1]/2+1)**2.0 ),
          x.shape[1]/2+1,axis=1), x.shape[0]/2+1,axis=0)+1e-9

    idt,dt=np.fft.ifft2,np.fft.fft2
    puOp=lambda x : idt( np.where(puRadius(x)==1e-9,1,puRadius(x)**-1.0)*dt(
          np.cos(x)*idt(puRadius(x)*dt(np.sin(x)))
         -np.sin(x)*idt(puRadius(x)*dt(np.cos(x))) ) )

    def phaseUnwrapper(ip):
       mirrored=np.zeros([x*2 for x in ip.shape])
       mirrored[:ip.shape[0],:ip.shape[1]]=ip
       mirrored[ip.shape[0]:,:ip.shape[1]]=ip[::-1,:]
       mirrored[ip.shape[0]:,ip.shape[1]:]=ip[::-1,::-1]
       mirrored[:ip.shape[0],ip.shape[1]:]=ip[:,::-1]

       return (ip+2*np.pi*
             np.round((puOp(mirrored).real[:ip.shape[0],:ip.shape[1]]-ip)
             /2/np.pi))

    i = 0
    if iterations < 1:
        interations = 1
    while i < iterations:
        i += 1
        phase = phaseUnwrapper(phase)

    return phase


if __name__ == '__main__':
    descriptiontext = "Pre-average data.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_file', help='Filename of input dataset')
    parser.add_argument('parmdb_file', help='Filename of input direction-independent selfcal instrument parmdb')
    parser.add_argument('input_colname', help='Name of input column to pre-average')
    parser.add_argument('output_colname', help='Name of output column')
    args = parser.parse_args()

    main(args.ms_file, args.parmdb_file, args.input_colname, args.output_colname)
