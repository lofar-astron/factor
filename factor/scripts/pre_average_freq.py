#! /usr/bin/env python
"""
Script to pre-average data using a sliding Gaussian kernel in frequency
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import sys
import os
from scipy.ndimage.filters import gaussian_filter1d as gfilter
import casacore.tables as pt


def main(ms_input, input_colname, output_data_colname, output_weights_colname,
    stddev_nchan=1):
    """
    Pre-average data using a sliding Gaussian kernel in frequency

    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile
    input_colname : str
        Name of the column in the MS from which the data are read
    output_data_colname : str
        Name of the column in the MS into which the averaged data are written
    output_weights_colname : str
        Name of the column in the MS into which the averaged data weights are written
    """
    ms = pt.table(ms_input, readonly=False, ack=False)
    ant1_list = ms.getcol('ANTENNA1')
    ant2_list = ms.getcol('ANTENNA2')
    data = ms.getcol(input_colname)
    weights = ms.getcol('WEIGHT_SPECTRUM')
    flags = ms.getcol('FLAG')

    flags[ np.isnan(data) ] = True # flag NaNs
    weights = weights * ~flags # set weight of flagged data to 0

    # Check that all NaNs are flagged
    if np.count_nonzero(np.isnan(data[~flags])) > 0:
        logging.error('NaNs in unflagged data in {0}!'.format(ms_input))
        sys.exit(1)

    #    Multiply every element of the data by the weights, convolve both
    #    the scaled data and the weights, and then divide the convolved data
    #    by the convolved weights (translating flagged data into weight=0).
    #    That's basically the equivalent of a running weighted average with
    #    a Gaussian window function.

    # weigth data and set bad data to 0 so nans do not propagate
    data = np.nan_to_num(data*weights)

    # smear weighted data and weights
    dataR = gfilter(np.real(data), stddev_nchan, axis=1)
    dataI = gfilter(np.imag(data), stddev_nchan, axis=1)
    weights = gfilter(weights, stddev_nchan, axis=1)

    # re-create data
    data = (dataR + 1j * dataI)
    data[(weights != 0)] /= weights[(weights != 0)] # avoid divbyzero

    # Add the output columns if needed
    if output_data_colname not in ms.colnames():
        desc = ms.getcoldesc(input_colname)
        desc['name'] = output_data_colname
        ms.addcols(desc)
    if output_weights_colname not in ms.colnames():
        desc = ms.getcoldesc('WEIGHT_SPECTRUM')
        desc['name'] = output_weights_colname
        ms.addcols(desc)

    ms.putcol(output_data_colname, data)
    ms.putcol('FLAG', flags) # this saves flags of nans, which is always good
    ms.putcol(output_weights_colname, weights)
    ms.close()
