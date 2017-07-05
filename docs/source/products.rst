.. _products:

Output
======

Directory structure
-------------------

Factor produces the following output inside the working directory:

``factor.log``
    Log file containing only the higher-level log messages. Detailed logs for each operation are available in the ``logs`` directory (see below).

``factor_directions.txt``
    Optional file listing the DDE calibrators. This file is generated only when no directions file is supplied by the user.

``chunks/``
    Directory containing the time-chunked datasets that Factor uses for processing.

``logs/``
    Directory containing the detailed operation logs.

    .. note::

        The log of each operation is stored as ``logs/operation_name/direction_name.out.log``. For example, the log of the ``facetselfcal`` operation for a direction named ``facet_patch_200`` will be stored in ``logs/facetselfcal/facet_patch_200.out.log``.

    .. note::

        Some error messages are stored in the ``logs/operation_name/direction_name.err.log`` file, but these are rarely of interest. Generally, important error messages will appear in the ``logs/operation_name/direction_name.out.log`` file. These log files can be very large, so a search for "error" is usually the easiest way to find any error messages.

``regions/``
    Directory containing the ds9 region files for the facet and self-calibration images. The following region files are made:

    * ``calimages_ds9.reg`` - the self-calibration image regions. The inner box shows the area over which sources are added back and cleaned. The outer box shows the area that is imaged.
    * ``facets_ds9.reg`` - the facet image regions.

``results/``
    Directory containing the results (images, etc.) of each operation. See :ref:`operations` for details of the primary output products of each operation.

    .. note::

        The output of each operation is stored in a directory named ``results/operation_name/direction_name/``. For example, the results of the ``facetselfcal`` operation for a direction named ``facet_patch_200`` will be stored in ``results/facetselfcal/facet_patch_200/``.

``state/``
    Directory containing files that save the state of a reduction.


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
