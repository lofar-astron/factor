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


def convert_fits_to_image(fitsimage):
    """
    Convert a fits image to a CASA image
    """
    import pyrap.images as pim

    outfilename = fitsimage.split('.fits')[0] + '.image'
    casaimage = pim.image(fitsimage)
    casaimage.saveas(outfilename, overwrite=True)

    return outfilename

