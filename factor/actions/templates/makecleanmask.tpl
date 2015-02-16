#!/usr/bin/env python

from lofar import bdsm
import pyrap.images
import numpy
import os
import pyfits
import sys

image_name = sys.argv[1]
mask_name  = image_name.split('.image')[0] + '.cleanmask'
atrous_do = {{ atrous_do}}
threshisl = {{ threshisl }}
threshpix = {{ threshpix }}

if atrous_do:
   threshisl = 4.0

img = bdsm.process_image(image_name, mean_map='zero', rms_box=(70,10),
    thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
    atrous_do=atrous_do, ini_method='curvature', quiet=True, stop_at='isl')

img.export_image(img_type='island_mask', mask_dilation=2, outfile=mask_name,
    img_format='casa', clobber=True)
