#! /usr/bin/env python
"""
Script to make a normal sky model from a polynomial sky model
"""
import argparse
from argparse import RawTextHelpFormatter
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import Angle
from astropy.table import Table
import numpy as np
from numpy.polynomial.polynomial import polyval
import casacore.tables as pt
import re
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


def RA2Angle(RA):
    """
    Returns Angle objects for input RA values.

    Parameters
    ----------
    RA : str, float or list of str, float
        Values of RA to convert. Can be strings in makesourcedb format or floats
        in degrees.

    Returns
    -------
    RAAngle : astropy.coordinates.Angle object

    """
    import astropy.units as u

    if type(RA) is not list:
        RA = [RA]

    if type(RA[0]) is str:
        try:
            RAAngle = Angle(Angle(RA, unit=u.hourangle), unit=u.deg)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise ValueError('RA not understood (must be string in '
                'makesourcedb format or float in degrees): {0}'.format(e.message))
    else:
        RAAngle = Angle(RA, unit=u.deg)

    return RAAngle


def Dec2Angle(Dec):
    """
    Returns Angle objects for input Dec values.

    Parameters
    ----------
    Dec : str, float or list of str, float
        Values of Dec to convert. Can be strings in makesourcedb format or floats
        in degrees

    Returns
    -------
    DecAngle : astropy.coordinates.Angle object

    """
    import astropy.units as u

    if type(Dec) is not list:
        Dec = [Dec]

    if type(Dec[0]) is str or type(Dec[0]) is np.string_:
        try:
            DecAngle = Angle(Dec, unit=u.deg)
        except KeyboardInterrupt:
            raise
        except ValueError:
            try:
                DecSex = [decstr.replace('.', ':', 2) for decstr in Dec]
                DecAngle = Angle(DecSex, unit=u.deg)
            except Exception as e:
                raise ValueError('Dec not understood (must be string in '
                    'makesourcedb format or float in degrees): {0}'.format(e.message))
        except Exception as e:
            raise ValueError('Dec not understood (must be string in '
                'makesourcedb format or float in degrees): {0}'.format(e.message))
    else:
        DecAngle = Angle(Dec, unit=u.deg)

    return DecAngle


def processLine(line, ncols):
    """
    Processes a makesourcedb line.

    Parameters
    ----------
    line : str
        Data line
    ncols : int
        Number of columns

    Returns
    -------
    line : str
        Processed line

    """
    if line.startswith("FORMAT") or line.startswith("format") or line.startswith("#"):
        return None

    # Check for SpectralIndex or SpectralTerms entries, which are unreadable as they use
    # the same separator for multiple orders as used for the columns
    line = line.strip('\n')
    a = re.search('\[.*\]', line)
    if a is not None:
        b = line[a.start(): a.end()]
        c = b.strip('[]')
        if ',' in c:
            c = c.replace(',', ';')
        line = line.replace(b, c)
    colLines = line.split(',')

    while len(colLines) < ncols:
        colLines.append(' ')

    return ','.join(colLines)


def processSpectralTerms(terms):
    if type(terms) is str or type(terms) is np.string_:
        terms = terms.strip('[]').split(';')
        terms = [float(t) for t in terms]
    return terms


def main(model_root, ms_file, skymodel, fits_mask=None, min_peak_flux_jy=0.0001,
    max_residual_jy=0.0):
    """
    Make a makesourcedb sky model for input MS from WSClean fits model images

    Parameters
    ----------
    model_root : str
        Root name of WSClean polynomial model sky model
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

    """
    min_peak_flux_jy = float(min_peak_flux_jy)
    max_residual_jy = float(max_residual_jy)

    if type(fits_mask) is str:
        if fits_mask.lower() == 'none':
            fits_mask = None

    # Read MS file and get the frequency info
    if '[' in ms_file and ']' in ms_file:
        files = ms_file.strip('[]').split(',')
        files = [f.strip() for f in files]
        ms_file = files[0]
    sw = pt.table(ms_file+'::SPECTRAL_WINDOW', ack=False)
    ms_freq = sw.col('REF_FREQUENCY')[0]
    sw.close()

    # Read in sky model
    ncols = 8
    outlines = ['Name, Type, Ra, Dec, SpectralTerms, MajorAxis, MinorAxis, Orientation']
    polymodel = model_root + '-components.txt'
    with open(polymodel) as f:
        for line in f:
            if line.startswith("# ReferenceFrequency"):
                ref_freq = float(line.split('=')[1].strip())
            else:
                outline = processLine(line, ncols)
                if outline is not None:
                    outlines.append(outline)
    data = Table.read(outlines, format='ascii', guess=False, delimiter=',',
        header_start=0, data_start=1)
    nsources = len(data)

    # Find model images and read in frequencies
    fits_models = glob.glob(model_root+'-00*-model.fits')
    if len(fits_models) == 0:
        # No channels images found, so look for non-MFS images
        fits_models = glob.glob(model_root+'-model.fits')
    if len(fits_models) == 0:
        print('ERROR: no model images found')
        sys.exit(1)
    freqs = []
    for f in fits_models:
        # Get the frequency info
        hdr = fits.getheader(f, 0, ignore_missing_end=True)
        freqs.append(hdr['CRVAL3']) # Hz

    # Determine nearest model frequency to MS frequency. We will use this to
    # get the fluxes (to match the frequency blocks used during imaging and
    # calibration)
    sky_freq = min(freqs, key=lambda x:abs(x-ms_freq))

    # Check if fits mask is empty
    if fits_mask is not None:
        if fits_mask.lower() == 'empty':
            # Handle case in which no sources were found during masking
            nsources = 0
            mask = None
        else:
            mask = fits.getdata(fits_mask, 0, ignore_missing_end=True)
            hdr = fits.getheader(fits_mask, 0, ignore_missing_end=True)
            w = wcs.WCS(hdr)
    else:
        mask = None

    # Discard components not in the mask, if given
    if mask is not None:
        pix = w.wcs_world2pix(np.array([RA2Angle(data['Ra'].tolist()),
            Dec2Angle(data['Dec'].tolist()), [0]*len(data), [0]*len(data)]).T, 0)
        not_in_mask = []
        for i, p in enumerate(pix):
            if mask[0, 0, int(round(p[1])), int(round(p[0]))] < 1:
                not_in_mask.append(i)
        data.remove_rows(not_in_mask)

    # Write sky model
    with open(skymodel, 'w') as outfile:
        outfile.write('FORMAT = Name, Type, Ra, Dec, I, Q, U, V, ReferenceFrequency, '
            'MajorAxis, MinorAxis, Orientation\n')
        for i, s in enumerate(data):
            polyterms = processSpectralTerms(s['SpectralTerms'])
            flux = polyval(sky_freq/ref_freq, polyterms)
            name = 'cc{}'.format(i)
            if s['Type'] == 'POINT':
                outfile.write('{0}, POINT, {1}, {2}, {3}, 0.0, 0.0, 0.0, {4}, , , \n'
                    .format(name, s['Ra'], s['Dec'], flux, ms_freq))
            elif s['Type'] == 'GAUSSIAN':
                outfile.write('{0}, GAUSSIAN, {1}, {2}, {3}, 0.0, 0.0, 0.0, {4}, {5}, {6}, {7}\n'
                    .format(name, s['Ra'], s['Dec'], flux, ms_freq,
                    s['MajorAxis'], s['MinorAxis'], s['Orientation']))


if __name__ == '__main__':
    descriptiontext = "Make a makesourcedb sky model from WSClean fits model images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fits_model_root', help='Root of model images')
    parser.add_argument('ms_file', help='Filename of MS for which sky model is to be made')
    parser.add_argument('skymodel', help='Filename of output sky model')
    parser.add_argument('-f' '--fits_mask', help='Filename of fits mask', type=str, default=None)
    parser.add_argument('-p', '--min_peak_flux_jy', help='Minimum absolute value of flux in Jy', type=float, default=0.0001)
    parser.add_argument('-r', '--max_residual_jy', help='Maximum acceptible total residual absolute flux in Jy', type=float, default=0.0)
    parser.add_argument('-i', '--interp', help='Interpolation method', type=str, default='linear')
    args = parser.parse_args()
    main(args.fits_model_root, args.ms_file, args.skymodel, fits_mask=args.fits_mask,
        min_peak_flux_jy=args.min_peak_flux_jy, max_residual_jy=args.max_residual_jy,
        interp=args.interp)
