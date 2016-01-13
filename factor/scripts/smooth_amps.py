#! /usr/bin/env python
"""
Script to smooth and normalize amplitude solutions
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import numpy
import os
import lofar.parmdb
import math
import scipy.signal
import shutil


def median_window_filter(ampl, half_window, threshold):

    ampl_tot_copy = numpy.copy(ampl)
    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]

    median_array  = scipy.signal.medfilt(sol,half_window*2-1)

    sol_flag = numpy.zeros(ndata+2*half_window, dtype=bool)
    sol_flag_val = numpy.zeros(ndata+2*half_window, dtype=bool)

    for i in range(half_window, half_window + ndata):
        # Compute median of the absolute distance to the median.
        window = sol[i-half_window:i+half_window+1]
        window_flag = sol_flag[i-half_window:i+half_window+1]
        window_masked = window[~window_flag]

        if len(window_masked) < math.sqrt(len(window)):
            # Not enough data to get accurate statistics.
            continue

        median = numpy.median(window_masked)
        q = 1.4826 * numpy.median(numpy.abs(window_masked - median))

        # Flag sample if it is more than 1.4826 * threshold * the
        # median distance away from the median.
        if abs(sol[i] - median) > (threshold * q):
            sol_flag[i] = True

        idx = numpy.where(sol == 0.0) # to remove 1.0 amplitudes
        sol[idx] = True

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
           ampl_tot_copy[i] = median_array[half_window+i] # fixed 2012
    return ampl_tot_copy


def main(instrument_name, instrument_name_smoothed):
    pol_list = ['0:0','1:1']
    gain = 'Gain'

    pdb = lofar.parmdb.parmdb(instrument_name)
    parms = pdb.getValuesGrid('*')

    key_names = parms.keys()
    nchans = len(parms[key_names[0]]['freqs'])

    # Get station names
    antenna_list = set([s.split(':')[-1] for s in pdb.getNames()])
    window = 4

    # Smooth
    for pol in pol_list:
        for antenna in antenna_list:
            for chan in range(nchans):
                real = numpy.copy(parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan])
                imag = numpy.copy(parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan])
                phase = numpy.arctan2(imag, real)
                amp = numpy.sqrt(imag**2 + real**2)

                amp = numpy.log10(amp)
                amp = median_window_filter(amp, window, 6)
                amp = median_window_filter(amp, window, 6)
                amp = median_window_filter(amp, 7, 6)
                amp = median_window_filter(amp, 4, 6)
                amp = median_window_filter(amp, 3, 6)
                amp = 10**amp

                # Clip extremely high amplitude solutions to prevent biasing the
                # normalization done below
                high_ind = numpy.where(amp > 5.0)
                amp[high_ind] = 5.0

                parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan] = amp * numpy.cos(phase)
                parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan] = amp * numpy.sin(phase)

    # Normalize the amplitude solutions to a mean of one across all channels
    amplist = []
    for chan in range(nchans):
        for pol in pol_list:
            for antenna in antenna_list:
                real = numpy.copy(parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan])
                imag = numpy.copy(parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan])
                amp  = numpy.copy(numpy.sqrt(real**2 + imag**2))
                amplist.append(amp)
    norm_factor = 1./(numpy.mean(amplist))

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
    args = parser.parse_args()

    main(args.instrument_name, args.instrument_name_smoothed)
