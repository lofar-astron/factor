#! /usr/bin/env python
"""
Script to reset amplitude solutions to unity
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import numpy
import os
import lofar.parmdb
import math
import shutil
import numpy as np


def main(instrument_name, instrument_name_reset):
    pdb = lofar.parmdb.parmdb(instrument_name)
    parms = pdb.getValuesGrid('*')
    key_names = parms.keys()

    # determine the number of polarizations in parmdb (2 or 4)
    if any('Gain:0:1:' in s for s in key_names):
        pol_list = ['0:0', '1:1', '0:1', '1:0']
    else:
        pol_list = ['0:0', '1:1']

    # Get station names
    antenna_list = list(set([s.split(':')[-1] for s in pdb.getNames()]))

    # Identify any gaps in time (frequency gaps are not allowed), as we need to handle
    # each section separately if gaps are present
    freqs = parms['Gain:1:1:Ampl:{s}'.format(s=antenna_list[0])]['freqs']
    freqwidths = parms['Gain:1:1:Ampl:{s}'.format(s=antenna_list[0])]['freqwidths']
    times = parms['Gain:1:1:Ampl:{s}'.format(s=antenna_list[0])]['times']
    timewidths = parms['Gain:1:1:Ampl:{s}'.format(s=antenna_list[0])]['timewidths']
    delta_times = times[1:] - times[:-1]
    gaps = np.where(delta_times > timewidths[:-1]*2.)
    gaps_ind = gaps[0] + 1
    gaps_ind = np.append(gaps_ind, np.array([len(times)]))

    # Reset the amplitude solutions to unity
    if os.path.exists(instrument_name_reset):
        shutil.rmtree(instrument_name_reset)
    pdbnew = lofar.parmdb.parmdb(instrument_name_reset, create=True)
    g_start = 0
    for g in gaps_ind:
        parms = pdb.getValues('*', freqs, freqwidths, times[g_start:g],
            timewidths[g_start:g], asStartEnd=False)
        for pol in pol_list:
            for antenna in antenna_list:
                try:
                    phase = np.copy(parms['Gain:'+pol+':Phase:'+antenna]['values'])
                    data_shape = phase.shape
                    pdbnew.addValues('Gain:'+pol+':Phase:{}'.format(antenna), phase, freqs, freqwidths,
                        times[g_start:g], timewidths[g_start:g], asStartEnd=False)
                    pdbnew.addValues('Gain:'+pol+':Ampl:{}'.format(antenna), numpy.ones(data_shape), freqs, freqwidths,
                        times[g_start:g], timewidths[g_start:g], asStartEnd=False)
                except KeyError:
                    # The antenna in question does not have solutions for this
                    # time period
                    continue
        g_start = g
    pdbnew.flush()


if __name__ == '__main__':
    descriptiontext = "Reset amplitude solutions to unity.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('instrument_name', help='name of the instrument parmdb to smooth')
    parser.add_argument('instrument_name_reset', help='name of the output parmdb')
    args = parser.parse_args()

    main(args.instrument_name, args.instrument_name_reset)
