from lofar import bdsm
import numpy

image_name = '{{ imagetomask }}'
atrous_do = {{ atrous_do}}
threshisl = {{ threshisl }}
threshpix = {{ threshpix }}
rmsbox = {{ rmsbox }}

if atrous_do:
   threshisl = 4.0

mask_name  = image_name.split('.image')[0] + '.cleanmask'

img = bdsm.process_image(image_name, mean_map='zero',
    thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
    atrous_do=atrous_do, ini_method='curvature', adaptive_rms_box=True,
    adaptive_thresh=20, quiet=True)

img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
    img_format='casa', clobber=True)
