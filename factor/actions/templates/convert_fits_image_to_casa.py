#!/usr/bin/env python
"""
Convert a fits image to a CASA image

For WSClean 1.7, use force_stokes_I = True to overide incorrect Stokes
keyword in model image headers. The restfreq is also set to avoid
problems with casapy2bbs.py
"""
import pyrap.images as pim
import pyrap.tables as pt
import numpy as np

if __name__ == '__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog <fits_filename> <casa_filename>')
    opt.add_option('-f', help='force output to be Stokes I', action='store_true',
        default=False)
    (options, args) = opt.parse_args()

    if len(args) != 2:
        opt.print_help()
        sys.exit(1)

    fitsimage = args[0]
    outfilename = args[1]
    casaimage = pim.image(fitsimage)
    casaimage.saveas(outfilename, overwrite=True)

    if option.f:
        coords = casaimage.coordinates().dict()
        coords['stokes1']['stokes'] = ['I']
        freq = coords['spectral2']['wcs']['crval']
        coords['spectral2']['restfreqs'] = np.array([freq])
        outtable = pt.table(outfilename, readonly=False, ack=False)
        outtable.putkeywords({'coords': coords})
        outtable.done()
