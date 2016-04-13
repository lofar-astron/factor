#! /usr/bin/env python
"""
Script to make a sky model for a facet
"""
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import aplpy
import numpy
import astropy.io.fits
import subprocess
import logging
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import NullFormatter


def meanclip(indata, clipsig=4.0, maxiter=10, converge_num=0.001, verbose=0):
   """
   Computes an iteratively sigma-clipped mean on a
   data set. Clipping is done about median, but mean
   is returned.
   .. note:: MYMEANCLIP routine from ACS library.
   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.
   Examples
   --------
   >>> mean, sigma = meanclip(indata)
   Parameters
   ----------
   indata: array_like
       Input data.
   clipsig: float
       Number of sigma at which to clip.
   maxiter: int
       Ceiling on number of clipping iterations.
   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.
   verbose: {0, 1}
       Print messages to screen?
   Returns
   -------
   mean: float
       N-sigma clipped mean.
   sigma: float
       Standard deviation of remaining pixels.
   """
   # Flatten array
   skpix = indata.reshape( indata.size, )

   ct = indata.size
   iter = 0; c1 = 1.0 ; c2 = 0.0

   while (c1 >= c2) and (iter < maxiter):
       lastct = ct
       medval = numpy.median(skpix)
       sig = numpy.std(skpix)
       wsm = numpy.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]

       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iter += 1
   # End of while loop

   mean  = numpy.mean( skpix )
   sigma = robust_sigma( skpix )

   if verbose:
       prf = 'MEANCLIP:'
       print '%s %.1f-sigma clipped mean' % (prf, clipsig)
       print '%s Mean computed in %i iterations' % (prf, iter)
       print '%s Mean = %.6f, sigma = %.6f' % (prf, mean, sigma)

   return mean, sigma


def find_imagenoise(data):
  mean, rms =  meanclip(data)
  return rms


def robust_sigma(in_y, zero=0):
    """
    Calculate a resistant estimate of the dispersion of
    a distribution. For an uncontaminated distribution,
    this is identical to the standard deviation.
    Use the median absolute deviation as the initial
    estimate, then weight points using Tukey Biweight.
    See, for example, Understanding Robust and
    Exploratory Data Analysis, by Hoaglin, Mosteller
    and Tukey, John Wiley and Sons, 1983.
    .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.
    Examples
    --------
    >>> result = robust_sigma(in_y, zero=1)
    Parameters
    ----------
    in_y : array_like
        Vector of quantity for which the dispersion is
        to be calculated
    zero : int
        If set, the dispersion is calculated w.r.t. 0.0
        rather than the central value of the vector. If
        Y is a vector of residuals, this should be set.
    Returns
    -------
    out_val : float
        Dispersion value. If failed, returns -1.
    """
    # Flatten array
    y = in_y.ravel()

    eps = 1.0E-20
    c1 = 0.6745
    c2 = 0.80
    c3 = 6.0
    c4 = 5.0
    c_err = -1.0
    min_points = 3

    if zero:
        y0 = 0.0
    else:
        y0 = numpy.median(y)

    dy    = y - y0
    del_y = abs( dy )

    # First, the median absolute deviation MAD about the median:

    mad = numpy.median( del_y ) / c1

    # If the MAD=0, try the MEAN absolute deviation:
    if mad < eps:
        mad = del_y.mean() / c2
    if mad < eps:
        return 0.0

    # Now the biweighted value:
    u  = dy / (c3 * mad)
    uu = u * u
    q  = numpy.where(uu <= 1.0)
    count = len(q[0])
    if count < min_points:
        module_logger.warn('ROBUST_SIGMA: This distribution is TOO WEIRD! '
                           'Returning {}'.format(c_err))
        return c_err

    numerator = numpy.sum( (y[q] - y0)**2.0 * (1.0 - uu[q])**4.0 )
    n    = y.size
    den1 = numpy.sum( (1.0 - uu[q]) * (1.0 - c4 * uu[q]) )
    siggma = n * numerator / ( den1 * (den1 - 1.0) )

    if siggma > 0:
        out_val = numpy.sqrt( siggma )
    else:
        out_val = 0.0

    return out_val


def main(imagefiles, maskfiles=None, imagenoise=None, interactive=False,
    facet_name=None):
    """
    Makes a png of the input images and masks

    Parameters
    ----------
    imagefiles : list or str
        Filenames of input image(s). These are assumed to follow the standard
        selfcal naming convention (i.e., with 'image02', 'image12', etc.)
    maskfiles : list or str, optional
        Filenames of input mask(s). If not given, these are assumed to follow
        the standard selfcal naming convention
    imagenoise : float, optional
        Image noise to use to set the image scaling
    interactive : bool, optional
        If True, plot images as grid and show. If False, save one image for
        each input image
    facet_name : str, optional
        Facet name for figure window

    """
     # Set logging level to ERROR to suppress extraneous info from aplpy
    logging.root.setLevel(logging.ERROR)

    if type(imagefiles) is str:
        imagefiles = imagefiles.strip('[]').split(',')
    imagefiles = [f.strip() for f in imagefiles]

    if maskfiles is None:
        maskfiles = []
        for imagefile in imagefiles:
            i = int(imagefile.split('casa_image')[1][0])
            maskfile = imagefile.split('.image')[0] + '.image.mask{}'.format(i)
            maskfile = maskfile.replace('image{}2'.format(i), 'image{}1'.format(i))
            if os.path.exists(maskfile):
                maskfiles.append(maskfile)
            else:
                maskfiles.append(None)
    else:
        if type(maskfiles) is str:
            maskfiles = maskfiles.strip('[]').split(',')
        maskfiles = [f.strip() for f in maskfiles]

    # Convert casa images to FITS if needed
    fitsfiles = []
    fitsmaskfiles = []
    for f, m in zip(imagefiles, maskfiles):
        if os.path.isdir(f):
            if not os.path.exists('{0}.fits'.format(f)):
                subprocess.call('image2fits in={0} out={0}.fits'.format(f),
                    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            fitsfiles.append('{}.fits'.format(f))
        else:
            fitsfiles.append(f)
        if m is not None:
            if os.path.isdir(m):
                if not os.path.exists('{0}.fits'.format(m)):
                    subprocess.call('image2fits in={0} out={0}.fits'.format(m),
                        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                fitsmaskfiles.append('{}.fits'.format(m))
            else:
                fitsmaskfiles.append(m)
        else:
            fitsmaskfiles.append(None)

    outplotname = fitsfiles[0].replace('.fits','.png')

    # find image noise
    if imagenoise is None:
        imagenoises = []
        for fitsimagename in fitsfiles:
            hdulist = astropy.io.fits.open(fitsimagename)
            data = hdulist[0].data
            imagenoises.append(find_imagenoise(data))
            hdulist.close()
        imagenoise = min(imagenoises)
    print('Image noise for plot scaling is: {} mJy'.format(imagenoise*1e3))

    # Set up plot(s)
    if interactive:
        row1_images = ['image02']
        row2_images = ['image12', 'image22']
        row3_images = ['image32', 'image42']
        num_tecamp_plots = 0
        for fitsimagename in fitsfiles:
            if any([fitsimagename.find(r) > 0 for r in row3_images]):
                num_tecamp_plots += 1
        Nc = 4
        Nr = 2 + int(numpy.ceil(num_tecamp_plots/4.0))
        fig = plt.figure(figsize=(14, Nr*3+1), facecolor='w', edgecolor='w')
        if facet_name is not None:
            fig.canvas.set_window_title('Selfcal Images for {0}'.format(facet_name))
        else:
            fig.canvas.set_window_title('Selfcal Images')
        gs = GridSpec(Nr, Nc, wspace=0.0, hspace=0.0)
        row1_colindx = 0
        row2_colindx = 0
        row3_colindx = 0
        first_gain = True

        for i, (fitsimagename, mask) in enumerate(zip(fitsfiles, fitsmaskfiles)):
            if any([fitsimagename.find(r) > 0 for r in row1_images]):
                ax = plt.subplot(gs[0, row1_colindx])
                if row1_colindx == 0:
                    plot_label = 'Dir. Indep.'
                else:
                    plot_label = None
                row1_colindx += 1
                subplotindx = row1_colindx
                for im in row1_images:
                    if im in fitsimagename:
                        if '_iter' in fitsimagename:
                            iter = int(fitsimagename.split('_iter')[1][0])
                            title = im + '_iter{}'.format(iter)
                        else:
                            title = im
            if any([fitsimagename.find(r) > 0 for r in row2_images]):
                ax = plt.subplot(gs[1, row2_colindx])
                if row2_colindx == 0:
                    plot_label = 'TEC'
                else:
                    plot_label = None
                row2_colindx += 1
                subplotindx = row2_colindx + 4
                for im in row2_images:
                    if im in fitsimagename:
                        if '_iter' in fitsimagename:
                            iter = int(fitsimagename.split('_iter')[1][0])
                            title = im + '_iter{}'.format(iter)
                        else:
                            title = im
            if any([fitsimagename.find(r) > 0 for r in row3_images]):
                if first_gain:
                    row_indx = 2
                if row3_colindx % 4 == 0:
                    row3_colindx = 0
                    if not first_gain:
                        row_indx += 1
                    else:
                        first_gain = False
                if row3_colindx == 0:
                    plot_label = 'TEC + Gain'
                else:
                    plot_label = None
                ax = plt.subplot(gs[row_indx, row3_colindx])
                row3_colindx += 1
                subplotindx = row3_colindx + 4 * row_indx
                for im in row3_images:
                    if im in fitsimagename:
                        if '_iter' in fitsimagename:
                            iter = int(fitsimagename.split('_iter')[1][0])
                            title = im + '_iter{}'.format(iter)
                        else:
                            title = im

            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')
            ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
            if plot_label is not None:
                ax.set_ylabel(plot_label)
            ax.xaxis.set_major_formatter(NullFormatter())
            ax.yaxis.set_major_formatter(NullFormatter())

            f = aplpy.FITSFigure(fitsimagename, figure=fig, slices=[0, 0],
                subplot=(Nr, Nc, subplotindx))
            f.show_colorscale(vmax=16*imagenoise, vmin=-6*imagenoise, cmap='bone')
            f.add_beam()
            f.tick_labels.hide()
            f.axis_labels.hide()
            f.set_title(title)
            f.beam.set_frame(True)
            f.beam.set_color('white')
            f.beam.set_edgecolor('black')
            f.beam.set_linewidth(1.)
            f.add_grid()
            f.grid.set_color('white')
            f.grid.set_alpha(0.5)
            f.grid.set_linewidth(0.2)
            if mask is not None:
                f.show_contour(mask, colors='red', levels=[0.0], filled=False, smooth=1,
                    alpha=0.6, linewidths=1)
        fig.show()
    else:
        for fitsimagename, mask in zip(fitsfiles, fitsmaskfiles):
            outplotname = fitsimagename.replace('.fits', '.png')
            f = aplpy.FITSFigure(fitsimagename, slices=[0, 0])
            f.show_colorscale(vmax=16*imagenoise, vmin=-6*imagenoise, cmap='bone')
            f.add_beam()
            f.beam.set_frame(True)
            f.beam.set_color('white')
            f.beam.set_edgecolor('black')
            f.beam.set_linewidth(1.)
            f.add_grid()
            f.grid.set_color('white')
            f.grid.set_alpha(0.5)
            f.grid.set_linewidth(0.2)
            f.add_colorbar()
            f.colorbar.set_axis_label_text('Flux (Jy beam$^{-1}$)')
            if mask is not None:
                f.show_contour(mask, colors='red', levels=[0.0], filled=False, smooth=1,
                    alpha=0.6, linewidths=1)
            f.save(outplotname, dpi=100, format='png')


if __name__ == '__main__':
    descriptiontext = "Make a sky model for a facet.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fullskymodel', help='name of the full skymodel')
    parser.add_argument('outmodel', help='name for the output')
    parser.add_argument('vertices_file', help='file containing facet vertices')
    parser.add_argument('-c', '--cal_only', help='return calibrator model only', type=bool, default=False)
    args = parser.parse_args()

    main(args.fullskymodel, args.outmodel, args.vertices_file, cal_only=args.cal_only)
