.. _products:

Primary output products
=======================

The primary output products of Factor are the images, solutions, sky models, and
calibrated data for one or more facets and the improved residual datasets of the
field. Each of these is described in detail below.


Final, primary-beam-corrected images
------------------------------------

The final images for the processed facets (which may consist only of the target
facet if ``image_target_only = True``) are mosaicked together and corrected for
the primary beam attenuation during the fieldmosaic operation (see
:ref:`fieldmosaic`). The fieldmosaic operation produces the following images (in
``results/fieldmosaic/field/``):

    * ``*correct_mosaic.pbcor.fits`` - the primary-beam-corrected mosaic. This image should be used for measurements of the source flux densities.
    * ``*correct_mosaic.pbcut.fits`` - the primary-beam-uncorrected (i.e., "flat-noise") mosaic (blanked beyond 40% power point of primary beam). This image can be used as the detection image for source detection in PyBDSF.

See :ref:`fieldmosaic` for more details.


Self-calibration solutions and sky models
-----------------------------------------

For each facet that undergoes self calibration, a set of self-calibration
solutions and sky models are generated. These solutions and models are used to
subtract or add sources in each facet. They are found in the facet directories
under ``results/facetselfcal``. See :ref:`facetselfcal` for more details.


Calibrated visibility data
--------------------------

For each facet that is imaged during the facetimage operation (see
:ref:`facetimage`), a set of measurement sets are made that contain the
calibrated visibilities for the facet. During the facetimage operation, the
sources for the given facet are first added back (sources in the field that lie
outside of the facet are not added back). Next, the visibilities are phase
shifted to the facet center, the calibration for the facet is applied, and
averaging is done (the amount of averaging allowed can be controlled with the
:term:`max_peak_smearing` parameter). These data are then used for the final
imaging of the facet.

If you want to image the facet outside of Factor (e.g., by hand with CASA or
WSClean), you can export these calibrated visibility data using archivefactor
with the ``-d`` flag. For example, to export the calibrated data for a facet
named ``facet_patch_468``, run the following::

    $ archivefactor factor.parset dir_output -d facet_patch_468

This command will create a directory inside ``dir_output`` named ``calibrated_data``,
inside of which will be measurement sets with the calibrated data (in the DATA
column). Note that these data have not been corrected for the primary-beam
attenuation but have had been corrected for the beam effects at the field phase
center.


Residual data
-------------

The subtraction of sources in the residual datasets (originally produced by
Prefactor) is improved during a Factor run. These datasets are those in the
``chunks/`` directory. Normally, these data are not needed once a Factor run for
a field is complete, but they are needed if the Factor run is later resumed. The
``archivefactor`` script will archive these datasets if the ``-r`` flag is used::

    $ archivefactor factor.parset dir_output -r
