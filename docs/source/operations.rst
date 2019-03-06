.. _operations:

Operations
==========

Most of the processing performed by Factor is done in "operations," which are sets of steps that are grouped together. Operations are generally defined depending on their ability to be run in parallel (e.g., multiple facetselfcal operations can be run in parallel but only one facetsub operation can be run at a time). The available operations and the primary data products of each are described in detail below.


.. _outlierpeel:

outlierpeel
-----------

This operation peels an outlier source and subtracts it. Since no self calibration is done during peeling, a sky model of the source must either be supplied by the user or must be one of those included in Factor (see ``factor/skymodels`` in the Factor installation for a list of available ones). To specify that a source is an outlier source and should be peeled, set :term:`outlier_source` to ``True`` in the directions file.

Primary products (in ``results/outlierpeel/direction_name/``):
    * ``combined_solutions_unnorm.h5`` - the (unnormalized) self-calibration solutions table used to subtract the source
    * ``*make_selfcal_plots*.png`` - plots of the self-calibration solutions


.. _facetpeel:

facetpeel
---------

This operation peels a facet calibrator. It is the same as the outlierpeel operation except that the facet is imaged (without adding back the calibrator) to improve the subtraction.

Primary products (in ``results/facetpeel/direction_name/``):
    * ``combined_solutions_unnorm.h5`` - the (unnormalized) self-calibration solutions table used to subtract the source
    * ``combined_solutions_norm.h5`` - the (normalized) self-calibration solutions table used to correct the facet data before imaging
    * ``*make_selfcal_plots*.png`` - plots of the self-calibration solutions
    * ``*image_full*`` - facet image made with part of the total bandwidth (not made if direction is a patch)
    * ``*wsclean_pre-image.fits`` - residual image of field for middle band before subtraction of new model
    * ``*wsclean_post-image.fits`` - residual image of field for middle band after subtraction of new model


.. _facetselfcal:

facetselfcal
------------

This operation self calibrates a facet calibrator and images the facet with part of the bandwidth.

.. note::

    The facet image is made only for facet-type directions (i.e., not the small patch-type directions that lie outside of the faceting radius). It is typically made with only a fraction of the total bandwidth (but distributed to sample the full bandwidth; see :term:`fractional_bandwidth_selfcal_facet_image` for details) and is used to improve the subtraction of non-calibrator sources in the facet. It is not the final facet image (which is made in the facetimage operation).

Primary products (in ``results/facetselfcal/direction_name/``):
    * ``combined_solutions.h5`` - the (normalized) self-calibration solutions table, with dTEC, CommonScalarPhase, and Gain solutions for all times and frequencies
    * ``*make_selfcal_plots*.png`` - plots of the self-calibration solutions
    * ``*wsclean_image*.fits`` - self-calibration images (FITS format)
    * ``*wsclean_image*.png`` - self-calibration images (png format)
    * ``*image_full*.fits`` - facet image made with part of the total bandwidth (not made if direction is a patch)
    * ``*wsclean_pre-image.fits`` - residual image of field for middle band before subtraction of new model
    * ``*wsclean_post-image.fits`` - residual image of field for middle band after subtraction of new model


.. _facetsub:

facetsub
--------

This operation subtracts the improved model of the facet or calibrator made in the outlierpeel, facetpeel, or facetselfcal operations.

Primary products:
    * None. The ``CORRECTED_DATA`` column of the chunked MS files in ``chunks/`` is updated.


.. _facetimage:

facetimage
----------

This operation images a facet using the full bandwidth. The self-calibration solutions from the same facet (if self calibration was successful) or the solutions from the nearest facet (if self calibration was unsuccessful) are used.

.. note::

    Multiple facetimage operations may be run if more than one set of imaging parameters were specified in the parset. In this case, the operation names will include the parameter values. E.g., ``facetimage_c15.0r-1.0t45.0u80.0`` is the operation to make an image with a cellsize of 15.0 arcsec, a robust value of -1.0, a Gaussian taper of 45.0 arcsec, and a uv cut of 80 lambda.

Primary products (in ``results/facetimage/direction_name/``):
    * ``*image_full*`` - facet image
    * files listed in ``mapfiles/imaging_input.mapfile`` - averaged, corrected uv data (kept only if :term:`keep_avg_facet_data` is ``True`` in the Factor parset)


.. _fieldmosaic:

fieldmosaic
-----------

This operation creates a mosaic of the field by joining the facet images together and correcting for the primary beam attenuation.

.. note::

    Multiple fieldmosaic operations may be run if more than one set of imaging parameters were specified in the parset. In this case, the operation names will include the parameter values. E.g., ``fieldmosaic_c15.0r-1.0t45.0u80.0`` is the operation to mosaic the images with a cellsize of 15.0 arcsec, a robust value of -1.0, a Gaussian taper of 45.0 arcsec, and a uv cut of 80 lambda.

Primary products (in ``results/fieldmosaic/field/``):
    * ``*correct_mosaic.pbcor.fits`` - the primary-beam-corrected mosaic. This image should be used for measurements of the source flux densities.
    * ``*correct_mosaic.pbcut.fits`` - the uncorrected mosaic (blanked beyond 40% power point of primary beam). This image can be used as the detection image for source detection in PyBDSF.


