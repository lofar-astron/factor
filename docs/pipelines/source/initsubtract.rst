.. _initial_subtract_operation:

Initial Subtract Operation
==========================

This section describes the initial subtract operation of Factor, denoted
``initsubtract``, which images each band at high and low resolution to make and subtract
sky models (two resolutions are used to capture the complete source structure and allow
imaging of the first side lobe of the primary beam).


.. _high_res_images:

Make high-res images
--------------------

Input
	MS files concatenated into bands of 10 subbands (2 MHz) each and their dir-independent parmdbs.

	The data should have clock offsets removed, the beam at the phase center
	applied (to ``DATA`` and ``CORRECTED_DATA`` columns), and be calibrated to
	25 arcsec resolution (in the ``CORRECTED_DATA`` column) using their dir-independent parmdbs.

Output
    For each band, a wide-field (~ 7 degree radius) image is made at
    approximately 25 arcsec resolution.

    .. note::

        The image should show the
        higher-resolution structure. There may be strong artifacts around bright
        sources.

    A region (~ 3 degrees across) of an example image is shown
    in Figure :num:`high-resolution-example-image`.

    .. _high-resolution-example-image:

    .. figure:: high_res_image.png
       :scale: 40 %
       :figwidth: 75 %
       :align: center
       :alt: example image

       High-resolution example image

Pipeline Steps
    create_ms_map, create_parmdb_map, create_high_sizes_map, create_low_sizes_map
        Make datamaps for input MS files, their dir-independent parmdbs, and
        the high- and low-res image sizes.

    wsclean_high1, mask_high, wsclean_high2
        High-res WSClean imaging run.


Make high-res sky models
------------------------

Input
    Model images created during high-res imaging (see :ref:`high_res_images`), one for each band.

Output
    Sky models in ``makesourcedb`` format, one for each band. A model from a typical 8-hour observation should
    contain ~ 10000 clean components.

Pipeline Steps
    fits_to_image_high
        Convert WSClean FITS model image to CASA model image.

    casa_to_bbs_high
        Convert CASA model image to ``makesourcedb`` sky model.


Subtract high-res models
------------------------

Input
    Model images created during high-res imaging (see :ref:`high_res_images`), one for each band.

Output
    ``SUBTRACTED_DATA`` and ``CORRECTED_SUBTRACTED_DATA`` columns for each band with all high-res sources subtracted.

Pipeline Steps
    create_model_high_map
        Make datamap for high-res model images.

    wsclean_ft_high
        Call WSClean to FT model image into MODEL_DATA column of each band.

    subtract_high
        Call BBS to subtract MODEL_DATA column from DATA column and correct subtracted column
        with dir-independent solutions.


.. _low_res_images:

Make low-res images
--------------------

Input
	Output of previous subtract step (``CORRECTED_SUBTRACTED_DATA`` columns)

Output
    For each band, a wide-field (~ 15 degree radius) image is made at
    approximately 75 arcsec resolution.

    .. note::

        The image should show the lower-resolution structure that was not
        picked up in the high-resolution images.

    A region (~ 3 degrees across) of an example image is shown
    in the `low-resolution-example-image`_.

    .. _low-resolution-example-image:

    .. figure:: low_res_image.png
       :scale: 40 %
       :figwidth: 75 %
       :align: center
       :alt: example image

       Low-resolution example image

Pipeline Steps
    average
        Average the ``CORRECTED_SUBTRACTED_DATA`` column as input to imager.

    wsclean_low1, mask_low, wsclean_low2
        Low-res WSClean imaging run.


Make low-res sky models
-----------------------

Input
    Model images created during low-res imaging (see :ref:`low_res_images`), one for each band.

Output
    Sky models in ``makesourcedb`` format, one for each band. A model from a typical 8-hour observation should
    contain ~ 1500 clean components.

Pipeline Steps
    fits_to_image_low
        Convert WSClean FITS model image to CASA model image.

    casa_to_bbs_low
        Convert CASA model image to ``makesourcedb`` sky model.


Subtract low-res models
------------------------

Input
    Model images created during low-res imaging (see :ref:`low_res_images`), one for each band.

Output
    ``SUBTRACTED_DATA_ALL`` column for each band with all low- and high-res sources subtracted.

Pipeline Steps
    create_model_low_map
        Make datamap for low-res model images.

    wsclean_ft_low
        Call WSClean to FT model image into MODEL_DATA column of each band.

    subtract_low
        Call BBS to subtract ``MODEL_DATA`` column from ``SUBTRACTED_DATA`` column.


Merge low- and high-res sky models
----------------------------------

Input
	Low- and high-res sky models in ``makesourcedb`` format, one of each for each band.

Output
    Merged sky models in ``makesourcedb`` format with both low- and high-res sources, one for each band.

Pipeline Steps
    merge
        Call LSMTool to merge low- and high-res sky models into a single sky model.


Partial initial subtract operation
----------------------------------

.. note::

    This step is done only for the ``initsubtract_subonly_pipeline.parset`` pipeline and replaces all of the above steps.

Input
    MS files concatenated into bands of 10 subbands (2 MHz) each, their dir-independent parmdbs, and their sky models.

Output
    ``SUBTRACTED_DATA_ALL`` column for each band with all low- and high-res sources subtracted.

Pipeline Steps
    create_ms_map, create_parmdb_map, create_skymodel_map
        Make datamaps for input MS files, their dir-independent parmdbs, and
        the sky models.

    subtract
        Call BBS to subtract the input sky models from the ``DATA`` column to make the ``SUBTRACTED_DATA`` column.



