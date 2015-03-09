from lofar import bdsm
import numpy
import sys
import os

ms        = sys.argv[6]
imageout  = sys.argv[7]
mask      = '{{ mask }}'
threshold = '{{ threshold }}'
threshold_5rms = threshold
use_rms = {{ use_rms }}
uvrange   = '{{ uvrange }}'
niter     = {{ niter }}
nterms    = {{ nterms }}
imsize    = [{{ imsize }}, {{ imsize }}]
mscale    = {{ mscale }}
cell      = ['{{ cell }}', '{{ cell }}']
cycfactor = {{ cycfactor }}
scales    = {{ scales }}
timer     = '{{ timer }}'
nfacets   = {{ nfacets }}
wplanes   = {{ wplanes }}
ncycles   = {{ ncycles }}
atrous_do = {{ atrous_do}}
threshisl = {{ threshisl }}
if atrous_do:
   threshisl = 4.0
threshpix = {{ threshpix }}
rmsbox = {{ rmsbox }}

# Change to imaging directory, since makemask cannot be used with absolute paths
dirname = os.path.dirname(imageout)
os.chdir(dirname)

for i in range(ncycles):
    # Set threshold of last cycle to 5 * rms
    if i == ncycles-1 and use_rms:
        threshold = threshold_5rms

    # Image for niter iterations
    clean(vis=ms,imagename=imageout,outlierfile="",field="",spw="",selectdata=True,timerange=timer,
          uvrange=uvrange,antenna="",scan="",observation="",mode="mfs",gridmode="widefield",wprojplanes=wplanes,
          facets=nfacets,cfcache="cfcache.dir",painc=360.0,epjtable="",interpolation="linear",
          niter=niter,gain=0.1,threshold=threshold,psfmode="clark",imagermode="csclean",
          ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=scales,negcomponent=-1,
          smallscalebias=0.6,interactive=False,mask=mask,nchan=-1,start=0,width=1,outframe="",
          veltype="radio",imsize=imsize,cell=cell,phasecenter="",restfreq="",stokes="I",
          weighting="briggs",robust=-0.25,uvtaper=False,outertaper=[''],innertaper=['1.0'],
          modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",
          npixels=0,npercycle=100,cyclefactor=cycfactor,cyclespeedup=-1,nterms=nterms,reffreq="",
          chaniter=False,flatnoise=True,allowchunk=False)

    # Make / refine clean mask and find noise. First export to FITS format, as
    # pyrap does not load in casapy
    if nterms > 1:
        image_name = imageout + '.image.tt0'
    else:
        image_name = imageout + '.image'
    fits_image = image_name + '.fits'
    exportfits(imagename=image_name, fitsimage=fits_image, overwrite=True)

    # Now run PyBDSM
    img = bdsm.process_image(fits_image, mean_map='zero',
        thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
        atrous_do=atrous_do, ini_method='curvature', adaptive_rms_box=True,
        adaptive_thresh=20, quiet=True)
    threshold_5rms = '{0}Jy'.format(img.clipped_rms*5.0)
    fits_mask  = imageout + '.cleanmask.fits'
    img.export_image(img_type='island_mask', mask_dilation=0, outfile=fits_mask,
        img_format='fits', clobber=True)

    # Now import FITS mask
    mask_image_temp = imageout + '.cleanmask.temp'
    importfits(fitsimage=fits_mask, imagename=mask_image_temp, overwrite=True)

    # Now match mask to image
    mask_image = imageout + '.cleanmask'
    makemask(mode='copy', inpimage=os.path.basename(image_name), inpmask=os.path.basename(mask_image_temp), output=os.path.basename(mask_image), overwrite=True)
    mask = mask_image

# Reimage from scratch with last mask
if image_final:
    mask_image_final = '
    os.system('cp -r {0} {1}'.format(mask_image, mask_image_final))
    os.system('rm -rf {0}*'.format(imageout))
    clean(vis=ms,imagename=imageout,outlierfile="",field="",spw="",selectdata=True,timerange=timer,
        uvrange=uvrange,antenna="",scan="",observation="",mode="mfs",gridmode="widefield",wprojplanes=wplanes,
        facets=nfacets,cfcache="cfcache.dir",painc=360.0,epjtable="",interpolation="linear",
        niter=niter,gain=0.1,threshold=threshold,psfmode="clark",imagermode="csclean",
        ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=scales,negcomponent=-1,
        smallscalebias=0.6,interactive=False,mask=mask,nchan=-1,start=0,width=1,outframe="",
        veltype="radio",imsize=imsize,cell=cell,phasecenter="",restfreq="",stokes="I",
        weighting="briggs",robust=-0.25,uvtaper=False,outertaper=[''],innertaper=['1.0'],
        modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",
        npixels=0,npercycle=100,cyclefactor=cycfactor,cyclespeedup=-1,nterms=nterms,reffreq="",
        chaniter=False,flatnoise=True,allowchunk=False)
