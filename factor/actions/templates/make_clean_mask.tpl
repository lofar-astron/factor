#!/usr/bin/env python

from lofar import bdsm
import pyrap.images
import numpy
import os
import pyfits

image_name = '{{ imagetomask }}'
atrous_do = {{ atrous_do}}
threshisl = {{ threshisl }}
threshpix = {{ threshpix }}
rmsbox = {{ rmsbox }}


if atrous_do:
   threshisl = 4.0
   print 'Changing island threshold to 4 because atrous_do=True'

mask_name  = image_name.split('.image')[0] + '.cleanmask'

print 'Making mask:', mask_name

img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
    thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
    atrous_do=atrous_do, ini_method='curvature', quiet=True)

img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
    img_format='casa', clobber=True)
