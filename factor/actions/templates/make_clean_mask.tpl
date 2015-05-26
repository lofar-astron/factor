from lofar import bdsm
import numpy
import sys

image_name = sys.argv[1]
mask_name  = sys.argv[2]
atrous_do = {{ atrous_do}}
threshisl = {{ threshisl }}
threshpix = {{ threshpix }}
rmsbox = {{ rmsbox }}
iterate_threshold = {{ iterate_threshold }}
adaptive_rmsbox = {{ adaptive_rmsbox }}
beam = {{ beam }}
img_format = '{{ format }}'

if atrous_do:
   threshisl = 4.0

if iterate_threshold:
    # Start with high threshold and lower it until we get at least one island
    nisl = 0
    threshpix = 25
    threshisl = 15
    while nisl == 0:
        img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
            thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
            atrous_do=atrous_do, ini_method='curvature', beam=beam,
            adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=20, quiet=True)
        nisl = img.nisl
        threshpix /= 1.2
        threshisl /= 1.2
else:
    img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
        thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
        atrous_do=atrous_do, ini_method='curvature', beam=beam,
        adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=20, quiet=True)

img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
    img_format=img_format, clobber=True)

log_file = mask_name + '.log'
with open(log_file, 'wb') as f:
    f.write('# 5-sigma clipped rms (Jy/beam): {0}'.format(5.0*img.clipped_rms))

