.. _facet_image:

Facet Final Imaging Operation
=============================

This operation performs imaging of the entire facet. The pipeline parset for this operation is ``facetimagefinal_pipeline.parset``.

.. note::

    There should be one pipeline per facet, and the pipelines may be run in parallel.

    This operation is separated from the add (:ref:`add_final_facet_sources`) operation as those pipelines
    must be run in series.


Data preparation
----------------

Input
	MS files from the :ref:`add_final_facet_sources` operation with phase-shifted facet sources in the
	``DATA`` column and the dir-dependent parmdbs from self calibration.

Output
    Datasets ready for imaging.

Pipeline Steps
    expand_merged_parmdb_map, apply_dir_dep
        Apply dir-dependent solutions from self calibration to the phase-shifted ``DATA`` column to make a ``CORRECTED_DATA`` column for imaging.

    average, create_compressed_mapfile, concat_averaged
        Average ``CORRECTED_DATA`` column in time and frequency and concatenate in frequency in preparation for imaging.

Test data
    With the phase-shifted facet all-source MS files (e.g., ``NEP_SB070-079.2ch10s.shift_all``) in ``Test_run/results/facetaddfinal/facet_patch_543/``, this step produces the MS file ``NEP_SB070-079.2ch10s.concat_averaged`` in ``Test_run/results/facetimagefinal/facet_patch_543/`` with averaged, concatenated (in frequency) ``DATA`` column.



Make image of entire facet
--------------------------

Input
	Full-resolution datasets (with all facet sources) with dir-dependent solutions applied.

Output
    Image of the entire facet. An example image is shown in Figure :num:`facet-example-image2`.

    .. note::

        The image should fully enclose the facet boundaries. Areas outside of the facet are not cleaned (and have all sources subtracted).

    .. _facet-example-image2:

    .. figure:: facet_image.png
       :scale: 80 %
       :figwidth: 75 %
       :align: center
       :alt: example image

       Facet example image

Pipeline Steps
    wsclean1, create_imagebase_map1, adjust_wsclean_mapfile1, copy_beam_info, mask, wsclean2, create_imagebase_map2
        WSClean imaging run. Imaging is done with a cell size of 1.5". Wide-band imaging is done if more than 5 bands are used. Multi-scale clean is not used, as WSClean does not currently support clean masks for this mode.

Test data
    With the the averaged, virtually-concatenated MS file ``RX42_SB070-079.2ch10s.concat_averaged``, this step produces the image ``NEP_SB070-079.2ch10s.wsclean2-image.fits`` (or ``NEP_SB070-079.2ch10s.wsclean2-MFS-image.fits`` if wide-band clean was used). All of these files are in ``Test_run/results/facetselfcal/facet_patch_543/``.


