#! /usr/bin/env python
"""
Script to make a sky model from fits model images
"""
import argparse
from argparse import RawTextHelpFormatter
from astropy.io import fits
from astropy import wcs
import numpy as np
import scipy.interpolate
import casacore.tables as pt
import sys
import os
import glob


def ra2hhmmss(deg):
    """Convert RA coordinate (in degrees) to HH MM SS"""

    from math import modf
    if deg < 0:
        deg += 360.0
    x, hh = modf(deg/15.)
    x, mm = modf(x*60)
    ss = x*60

    return (int(hh), int(mm), ss)


def dec2ddmmss(deg):
    """Convert DEC coordinate (in degrees) to DD MM SS"""

    from math import modf
    sign = (-1 if deg < 0 else 1)
    x, dd = modf(abs(deg))
    x, ma = modf(x*60)
    sa = x*60

    return (int(dd), int(ma), sa, sign)


def convert_radec_str(ra, dec):
    """Takes ra, dec in degrees and returns makesourcedb strings"""
    ra = ra2hhmmss(ra)
    sra = str(ra[0]).zfill(2)+':'+str(ra[1]).zfill(2)+':'+str("%.3f" % (ra[2])).zfill(6)
    dec = dec2ddmmss(dec)
    decsign = ('-' if dec[3] < 0 else '+')
    sdec = decsign+str(dec[0]).zfill(2)+'.'+str(dec[1]).zfill(2)+'.'+str("%.3f" % (dec[2])).zfill(6)

    return sra, sdec


def main(fits_model_root, ms_file, skymodel, fits_mask=None, min_peak_flux_jy=0.0001,
    max_residual_jy=0.05, interp='linear'):
    """
    Make a makesourcedb sky model for input MS from WSClean fits model images

    Parameters
    ----------
    fits_model_root : str
        Root name of WSClean fits model files (without the "-XXXX-model.fits" part)
    ms_file : str
        Filename of MS for which sky model is to be made. Can be a list of files
        (e.g., '[ms1,ms2,...]', in which case they should all have the same
        frequency
    skymodel : str
        Filename of the output makesourcedb sky model
    fits_mask : str, optional
        Filename of fits mask
    min_peak_flux_jy : float, optional
        Minimum absolute value of flux in Jy of a source in lowest-frequency model image
        to include in output model
    max_residual_jy : float, optional
        Maximum acceptible total residual absolute flux in Jy
    interp : str, optional
        Interpolation method. Can be any supported by scipy.interpolate.interp1d:
            'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'

    """
    min_peak_flux_jy = float(min_peak_flux_jy)
    max_residual_jy = float(max_residual_jy)

    if type(fits_mask) is str:
        if fits_mask.lower() == 'none':
            fits_mask = None

    # Find model images: look first for channel images and MFS image
    fits_models = glob.glob(fits_model_root+'-00*-model.fits')
    if len(fits_models) > 0:
        # Get the MFS image
        mfs_model = fits_model_root+'-MFS-model.fits'
    else:
        # No channels images found, so look for non-MFS images
        fits_models = glob.glob(fits_model_root+'-model.fits')
        mfs_model = None
    if len(fits_models) == 0:
        print('ERROR: no model images found')
        sys.exit(1)

    # Read MS file and get the frequency info
    if '[' in ms_file and ']' in ms_file:
        files = ms_file.strip('[]').split(',')
        files = [f.strip() for f in files]
        ms_file = files[0]
    sw = pt.table(ms_file+'::SPECTRAL_WINDOW', ack=False)
    ms_freq = sw.col('REF_FREQUENCY')[0]
    sw.close()

    # Read in model images
    freqs = []
    model_images = []
    for f in fits_models:
        # Get the frequency info
        hdr = fits.getheader(f, 0, ignore_missing_end=True)
        freqs.append(hdr['CRVAL3']) # Hz
        model_images.append(fits.getdata(f, 0, ignore_missing_end=True))
    w = wcs.WCS(hdr)

    # Read in MFS image
    if mfs_model is None:
        mfs_model = fits_models[0]
    mfs_image = fits.getdata(mfs_model, 0, ignore_missing_end=True)

    # Sort by freq
    sorted_ind = np.argsort(freqs)
    freqs = np.array(freqs)[sorted_ind]
    fits_models = np.array(fits_models)[sorted_ind]
    model_images = np.array(model_images)[sorted_ind]

    # Find pixels that meet the flux cut (and are in the mask, if given)
    if fits_mask is not None:
        if fits_mask.lower() == 'empty':
            # Handle case in which no sources were found during masking
            nonzero_ind = [[], []]
        else:
            mask = fits.getdata(fits_mask, 0, ignore_missing_end=True)
            nonzero_ind = np.where((np.abs(mfs_image) > min_peak_flux_jy) & (mask > 0))
    else:
        nonzero_ind = np.where(np.abs(mfs_image) > min_peak_flux_jy)

    # Interpolate the fluxes to the frequency of the MS
    nsources = len(nonzero_ind[0])
    fluxes = []
    names = []
    ras = []
    decs = []
    for i in range(nsources):
        index = [nonzero_ind[j][i] for j in range(4)]
        index.reverse() # change to WCS coords
        ras.append(w.wcs_pix2world(np.array([index]), 0, ra_dec_order=True)[0][0])
        decs.append(w.wcs_pix2world(np.array([index]), 0, ra_dec_order=True)[0][1])
        names.append('cc{}'.format(i))
        index.reverse() # change back to image coords
        flux_array = np.array([im[tuple(index)] for im in model_images])

        # If MS frequency lies outside range, just use nearest freq
        if ms_freq < freqs[0]:
            flux = flux_array[0]
        elif ms_freq > freqs[-1]:
            flux = flux_array[-1]
        else:
            # Otherwise interpolate
            flux = scipy.interpolate.interp1d(freqs, flux_array, kind=interp)(ms_freq)
        fluxes.append(flux)

    # Remove sources until we reach the desired residual
    if len(fluxes) > 0:
        total_flux = np.sum(np.abs(fluxes))
        keep_ind = np.where(np.abs(fluxes) > min_peak_flux_jy)
        while (total_flux - np.sum(np.abs(np.array(fluxes)[keep_ind]))) < max_residual_jy:
            min_peak_flux_jy *= 1.1
            keep_ind = np.where(np.abs(fluxes) > min_peak_flux_jy)
            if len(keep_ind[0]) < 50:
                # keep up to 50 sources regardless of the residual
                break
        fluxes = np.array(fluxes)[keep_ind]
        ras = np.array(ras)[keep_ind]
        decs = np.array(decs)[keep_ind]
        names = np.array(names)[keep_ind]

    # Write sky model
    with open(skymodel, 'w') as outfile:
        outfile.write('FORMAT = Name, Type, Ra, Dec, I, Q, U, V, ReferenceFrequency\n')
        for name, ra, dec, flux in zip(names, ras, decs, fluxes):
            ra_str, dec_str = convert_radec_str(ra, dec)
            outfile.write('{0}, POINT, {1}, {2}, {3}, 0.0, 0.0, 0.0, {4}\n'
                .format(name, ra_str, dec_str, flux, ms_freq))


if __name__ == '__main__':
    descriptiontext = "Make a makesourcedb sky model from WSClean fits model images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fits_model_root', help='Root of model images')
    parser.add_argument('skymodel', help='Filename of output sky model')
    args = parser.parse_args()
    main(args.fits_model_root, args.skymodel)
