from lofar import bdsm
import numpy

image_name = sys.argv[1]
mask_name  = image_name.split('.image')[0] + '.cleanmask'
atrous_do = {{ atrous_do}}
threshisl = {{ threshisl }}
threshpix = {{ threshpix }}
rmsbox = {{ rmsbox }}

if atrous_do:
   threshisl = 4.0

img = bdsm.process_image(image_name, mean_map='zero',
    thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
    atrous_do=atrous_do, ini_method='curvature', adaptive_rms_box=True,
    adaptive_thresh=20, quiet=True)

img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
    img_format='casa', clobber=True)

log_file = mask_name + '.log'
with open(log_file, 'wb') as f:
    f.write('# Clipped rms (Jy/beam): {0}'.format(img.clipped_rms))

