"""
Some functions called by multiple actions
"""


def make_basename(prefix, direction=None, band=None, index=None):
    """
    Returns a standard name pattern

    Parameters
    ----------
    prefix : str
        A prefix for the name, usually related to the particular operation step
    direction : Direction object, optional
        A direction
    band : Band object, optional
        A band
    index : int, optional
        An index for the particular operation step
    """
    if direction is not None:
        dirtxt = '_{0}'.format(direction.name)
    else:
        dirtxt = ''
    if band is not None:
        bandtxt = '_{0}'.format(band.name)
    else:
        bandtxt = ''
    if index is not None:
        indtxt = '-{0}'.format(index)
    else:
        indtxt = ''

    return '{0}{1}{2}{3}'.format(prefix, dirtxt, bandtxt, indtxt)


def make_image_basename(input_datamap, direction=None, band=None, prefix=None):
    """
    Define a standard name pattern for imaging files

    Parameters
    ----------
    input_datamap : Datamap
        Data map to files
    direction : Direction object, optional
        A direction
    band : Band object, optional
        A band
    prefix : str, optional
        String to prepend to the image basename. If None, 'image' is used

    """
    from factor.lib.datamap_lib import read_mapfile
    import re
    import os

    if prefix is None:
        prefix = ''
    else:
        prefix += '-'

    msfiles, hosts = read_mapfile(input_datamap)
    image_basenames = []

    for msfile in msfiles:
        msbase = os.path.basename(msfile)

        if direction is not None:
            dirtxt = '_{0}'.format(direction.name)
        else:
            dirtxt = ''
        if band is not None:
            bandtxt = '_{0}'.format(band.name)
        else:
            bandtxt = ''

        image_basenames.append('%s%s%s%s' % (prefix, re.sub(r'.MS|.ms', '', msbase), dirtxt, bandtxt))

    return image_basenames


def copy_model_images(modelbasename, ms_list, nterms=1):
    """
    Copies model images and returns list of output basenames
    """
    inmodelimages = []
    outmodelimages = []
    outmodelbasenames = []

    for i, band in enumerate(bands):
        outmodelbasenames.append(modelbasename + '_band{0}'.format(i))
        if nterms == 1:
            inmodelimages.append([modelbasename + '.model'])
            outmodelimages.append([modelbasename + '_band{0}.model'.format(i)])
        elif nterms == 2:
            inmodelimages.append([modelbasename + '.model.tt0',
                modelbasename + '.model.tt1'])
            outmodelimages.append([modelbasename + '_band{0}.model.tt0'.format(i),
                modelbasename + '_band{0}.model.tt1'.format(i)])
        else:
            inmodelimages.append([modelbasename + '.model.tt0',
                modelbasename + '.model.tt1', modelbasename + '.model.tt2'])
            outmodelimages.append([modelbasename + '_band{0}.model.tt0'.format(i),
                modelbasename + '_band{0}.model.tt1'.format(i),
                modelbasename + '_band{0}.model.tt2'.format(i)])

    for inmods, outmods in zip(inmodelimages, outmodelimages):
        for inmod, outmod in zip(inmods, outmods):
            if os.path.exists(outmod):
                os.system('rm -rf {0}'.format(outmod))
            os.system('cp -r {0} {1}'.format(inmod, outmod))

    return outmodelbasenames


def getOptimumSize(size):
    """
    Gets the nearest optimum image size

    Taken from the casa source code (cleanhelper.py)

    Parameters
    ----------
    size : int
        Target image size in pixels

    Returns
    -------
    optimum_size : int
        Optimum image size nearest to target size

    """
    import numpy

    def prime_factors(n, douniq=True):
        """ Return the prime factors of the given number. """
        factors = []
        lastresult = n
        sqlast=int(numpy.sqrt(n))+1
        if n == 1:
            return [1]
        c=2
        while 1:
             if (lastresult == 1) or (c > sqlast):
                 break
             sqlast=int(numpy.sqrt(lastresult))+1
             while 1:
                 if(c > sqlast):
                     c=lastresult
                     break
                 if lastresult % c == 0:
                     break
                 c += 1

             factors.append(c)
             lastresult /= c

        if (factors==[]): factors=[n]
        return  numpy.unique(factors).tolist() if douniq else factors

    n = int(size)
    if (n%2 != 0):
        n+=1
    fac=prime_factors(n, False)
    for k in range(len(fac)):
        if (fac[k] > 7):
            val=fac[k]
            while (numpy.max(prime_factors(val)) > 7):
                val +=1
            fac[k]=val
    newlarge=numpy.product(fac)
    for k in range(n, newlarge, 2):
        if ((numpy.max(prime_factors(k)) < 8)):
            return k
    return newlarge


def convert_fits_to_image(fitsimage, force_stokes_I=True):
    """
    Convert a fits image to a CASA image

    For WSClean 1.7, use force_stokes_I = True to overide incorrect Stokes
    keyword in model image headers. The restfreq is also set to avoid
    problems with casapy2bbs.py
    """
    import pyrap.images as pim
    import pyrap.tables as pt
    import numpy as np

    outfilename = fitsimage.split('.fits')[0] + '.image'
    casaimage = pim.image(fitsimage)
    casaimage.saveas(outfilename, overwrite=True)

    if force_stokes_I:
        coords = casaimage.coordinates().dict()
        coords['stokes1']['stokes'] = ['I']
        freq = coords['spectral2']['wcs']['crval']
        coords['spectral2']['restfreqs'] = np.array([freq])
        outtable = pt.table(outfilename, readonly=False, ack=False)
        outtable.putkeywords({'coords': coords})
        outtable.done()

    return outfilename


def get_val_from_str(val_str, retunits):
    """
    Convert a string with optional units to value in units of retunits

    If a string does not have any units, it's assumed to be in units of retunits
    """
    from itertools import groupby
    from astropy import units as u

    parts = [''.join(g).strip() for _, g in groupby(val_str, str.isalpha)]
    val = float(parts[0])
    if len(parts) > 1:
        if type(parts[1]) is str:
            if parts[1].lower() == 'e':
                # Check if number uses exponential notation (e.g., 1e8)
                parts = [parts[0] + parts[1] + parts[2]] + parts[3:]
            units = parts[1]
        else:
            units = None
    if units is None:
        units = retunits

    q = u.Quantity(val, units).to(retunits)

    return q.value
