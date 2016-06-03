.. _operations:

Operations
==========

Most of the processing performed by Factor is done in "operations," which are sets of steps that are grouped together. Operations are generally defined depending on their ability to be run in parallel (e.g., multiple facetselfcal operations can be run in parallel but only one facetsub operation can be run at a time). The available operations are described in detail below.


outlierpeel
-----------

This operation peels an outlier source and subtracts it. During peeling, no self calibration is done and hence a sky model of the source must either be supplied by the user or must be one of those included in Factor (see ``factor/skymodels`` in the Factor installation for a list of available ones). To specify that a source is an outlier source and should be peeled, set ``outlier_source = True`` in the directions file.

Primary products (in ``results/outlierpeel/direction_name/``):
    * ``*merge_selfcal_parmdbs`` - the (unnormalized) self-calibration solutions table used to subtract the source
    * ``*make_selfcal_plots*.png`` - plots of the self-calibration solutions


facetpeel
---------

This operation peels a facet calibrator. It is the same as the outlierpeel operation except that the facet is later imaged (without adding back the calibrator).

Primary products (in ``results/facetpeel/direction_name/``):
    * ``*merge_selfcal_parmdbs`` - the (unnormalized) self-calibration solutions table used to subtract the source
    * ``*merge_normalized_selfcal_parmdbs`` - the (normalized) self-calibration solutions table used to correct the facet data before imaging
    * ``*make_selfcal_plots*.png`` - plots of the self-calibration solutions (see Figure XX for an example)


facetselfcal
------------

This operation self calibrates a facet calibrator and images the facet.

.. note::

    The facet image is made only for facet-type directions (i.e., not patch-type directions).

Primary products (in ``results/facetselfcal/direction_name/``):
    * ``*merge_selfcal_parmdbs`` - the (normalized) self-calibration solutions table
    * ``*make_selfcal_plots*.png`` - plots of the self-calibration solutions
    * ``*casa_image*.image(.tt0)`` - self-calibration images (CASA format)
    * ``*casa_image*.png`` - self-calibration images (png format)
    * ``*image_full2*`` - facet image (not made if direction is a patch)
    * ``*wsclean_pre-image.fits`` - residual image of field for middle band before subtraction of new model
    * ``*wsclean_post-image.fits`` - residual image of field for middle band after subtraction of new model
    * files listed in ``mapfiles/concat_averaged_compressed.mapfile`` - averaged, corrected uv data (kept only if ``keep_avg_facet_data = True`` in the Factor parset)


facetsub
--------

This operation subtracts the improved model of the facet or calibrator made in the facetselfcal operation.

Primary products:
    * None. The ``SUBTRACTED_DATA_ALL_NEW`` column of the chunked MS files in ``chunks/`` are updated.


facetimage
----------

This operation images a facet using either the self-calibration solutions from the same facet (if self calibration was successful) or the solutions from a nearby facet (if self calibration was unsuccessful).

Primary products (in ``results/facetimage/direction_name/``):
    * ``*image_full2*`` - facet image
    * files listed in ``mapfiles/concat_averaged_compressed.mapfile`` - averaged, corrected uv data (kept only if ``keep_avg_facet_data = True`` in the Factor parset)


facetpeelimage
--------------

This operation images a facet using the self-calibration solutions from the facetpeel operation.

Primary products (in ``results/facetpeelimage/direction_name/``):
    * ``*image_full2*`` - facet image
    * files listed in ``mapfiles/concat_averaged_compressed.mapfile`` - averaged, corrected uv data (kept only if ``keep_avg_facet_data = True`` in the Factor parset)


fieldmosaic
-----------

This operation creates a mosaic of the field by joining the facet images together and correcting for the primary beam attenuation.

Primary products (in ``results/fieldmosaic/field/``):
    * ``*correct_mosaic.pbcor.fits`` - primary-beam-corrected mosaic
    * ``*correct_mosaic.pbcut.fits`` - primary-beam-corrected mosaic (blanked beyond 40% power point of primary beam)


