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


def main(instrument_name, instrument_name_reset):
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

    # Get station names
    antenna_list = set([s.split(':')[-1] for s in pdb.getNames()])

    # Reset the amplitude solutions to unity
    for chan in range(nchans):
        for pol in pol_list:
            for antenna in antenna_list:
                real = numpy.copy(parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan])
                imag = numpy.copy(parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan])
                phase = numpy.arctan2(imag, real)
                parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan] = numpy.copy(numpy.cos(phase))
                parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan] = numpy.copy(numpy.sin(phase))

    if os.path.exists(instrument_name_reset):
        shutil.rmtree(instrument_name_reset)
    pdbnew = lofar.parmdb.parmdb(instrument_name_reset, create=True)
    pdbnew.addValues(parms)
    pdbnew.flush()


if __name__ == '__main__':
    descriptiontext = "Reset amplitude solutions to unity.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('instrument_name', help='name of the instrument parmdb to smooth')
    parser.add_argument('instrument_name_reset', help='name of the output parmdb')
    args = parser.parse_args()

    main(args.instrument_name, args.instrument_name_reset)
