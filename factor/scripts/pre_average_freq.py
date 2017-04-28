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
from scipy.special import erf
import casacore.tables as pt
import pickle
import itertools


def main(ms_input, input_colname, output_data_colname, output_weights_colname,
    baseline_file, delta_theta_deg, target_peak_reduction_factor=0.99):
    """
    Pre-average data using a sliding Gaussian kernel in frequency

    Parameters
    ----------
    ms_input : str
        MS filename
    input_colname : str
        Name of the column in the MS from which the data are read
    output_data_colname : str
        Name of the column in the MS into which the averaged data are written
    output_weights_colname : str
        Name of the column in the MS into which the averaged data weights are
        written
    baseline_file : str
        Filename of pickled baseline lengths
    delta_theta_deg : float
        Radius of calibration region in degrees
    target_peak_reduction_factor : float, optional
        Target reduction in peak flux density. Note: this reduction is in
        addition to any incurred by earlier averaging

    """
    if os.path.exists(baseline_file):
        f = open(baseline_file, 'r')
        baseline_dict = pickle.load(f)
        f.close()
    else:
        print('Cannot find baseline_file. Exiting...')
        sys.exit(1)
    delta_theta_deg = float(delta_theta_deg)
    target_peak_reduction_factor = float(target_peak_reduction_factor)

    ms = pt.table(ms_input, readonly=False, ack=False)
    ant1_list = ms.getcol('ANTENNA1')
    ant2_list = ms.getcol('ANTENNA2')
    data_all = ms.getcol(input_colname)
    weights_all = ms.getcol('WEIGHT_SPECTRUM')
    flags = ms.getcol('FLAG')

    # Get lowest frequency of MS and channel width
    sw = pt.table(ms_input+'::SPECTRAL_WINDOW', ack=False)
    freq_hz = sw.col('CHAN_FREQ')[0][0]
    chan_width_hz = sw.col('CHAN_WIDTH')[0][0]

    flags[ np.isnan(data_all) ] = True # flag NaNs
    weights_all = weights_all * ~flags # set weight of flagged data to 0

    # Check that all NaNs are flagged
    if np.count_nonzero(np.isnan(data_all[~flags])) > 0:
        logging.error('NaNs in unflagged data in {0}!'.format(ms_input))
        sys.exit(1)

    # Weight data and set bad data to 0 so nans do not propagate
    data_all = np.nan_to_num(data_all*weights_all)

    # Iteration on baseline combination
    for ant in itertools.product(set(ant1_list), set(ant2_list)):
        if ant[0] >= ant[1]:
            continue
        sel1 = np.where(ant1_list == ant[0])[0]
        sel2 = np.where(ant2_list == ant[1])[0]
        sel_list = sorted(list(frozenset(sel1).intersection(sel2)))

        data = data_all[sel_list,:,:]
        weights = weights_all[sel_list,:,:]

        # compute the Gaussian sigma from the max bandwidth over which we
        # can average and avoid significant bandwidth smearing but limited to
        # no more than 3 MHz (to avoid smoothing over the beam-induced effects)
        lambda_km = 299792.458 / freq_hz
        dist_km = baseline_dict['{0}-{1}'.format(ant[0], ant[1])]
        resolution_deg = lambda_km / dist_km * 180.0 / np.pi
        stddev_hz = min(3e6, get_target_bandwidth(freq_hz, delta_theta_deg,
            resolution_deg, target_peak_reduction_factor)/4.0)
        stddev_nchan = stddev_hz / chan_width_hz * np.sqrt(0.5 / dist_km)

        # smear weighted data and weights
        dataR = gfilter(np.real(data), stddev_nchan, axis=1)
        dataI = gfilter(np.imag(data), stddev_nchan, axis=1)
        weights = gfilter(weights, stddev_nchan, axis=1)

        # re-create data
        data = (dataR + 1j * dataI)
        data[(weights != 0)] /= weights[(weights != 0)] # avoid divbyzero
        data_all[sel_list,:,:] = data
        weights_all[sel_list,:,:] = weights

    # Add the output columns if needed
    if output_data_colname not in ms.colnames():
        desc = ms.getcoldesc(input_colname)
        desc['name'] = output_data_colname
        ms.addcols(desc)
    if output_weights_colname not in ms.colnames():
        desc = ms.getcoldesc('WEIGHT_SPECTRUM')
        desc['name'] = output_weights_colname
        ms.addcols(desc)

    ms.putcol(output_data_colname, data_all)
    ms.putcol('FLAG', flags) # this saves flags of nans, which is always good
    ms.putcol(output_weights_colname, weights_all)
    ms.close()


def get_bandwidth_smearing_factor(freq, delta_freq, delta_theta, resolution):
    """
    Returns peak flux density reduction factor due to bandwidth smearing

    Parameters
    ----------
    freq : float
        Frequency at which averaging will be done
    delta_freq : float
        Bandwidth over which averaging will be done
    delta_theta : float
        Distance from phase center
    resolution : float
        Resolution of restoring beam

    Returns
    -------
    reduction_factor : float
        Ratio of post-to-pre averaging peak flux density

    """
    beta = (delta_freq/freq) * (delta_theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    reduction_factor = ((np.pi**0.5)/(gamma * beta)) * (erf(beta*gamma/2.0))

    return reduction_factor


def get_target_bandwidth(freq, delta_theta, resolution, reduction_factor):
    """
    Returns the bandwidth for given peak flux density reduction factor

    Parameters
    ----------
    freq : float
        Frequency at which averaging will be done
    delta_theta : float
        Distance from phase center
    resolution : float
        Resolution of restoring beam
    reduction_factor : float
        Ratio of post-to-pre averaging peak flux density

    Returns
    -------
    delta_freq : float
        Bandwidth over which averaging will be done
    """
    # Increase delta_freq until we drop below target reduction_factor
    delta_freq = 1e-3 * freq
    while get_bandwidth_smearing_factor(freq, delta_freq, delta_theta,
        resolution) > reduction_factor:
        delta_freq *= 1.1

    return delta_freq
