.. _facet_image:

Facet Imaging Operation
=============================

This section describes the facet imaging operation of Factor, denoted
``facetimage``, which performs imaging of the entire facet when either the facet did not successfully go through selfcal or reimaging is desired. There are two possible pipeline parsets for this operation:

Reimaging pipeline
    For use when the facet under consideration has successfully gone through self calibration and reimaging is desired.

Non-reimaging pipeline
    For use when the facet under consideration has not successfully gone through self calibration. In this case, the self-calibration solutions from the nearset facet should be used instead.

.. note::

    There should be one pipeline per facet, and the pipelines may be run in parallel.


Add sources from CC model
-------------------------

.. note::

    This step is done only for non-reimaging pipeline.

Input
	MS files from the :ref:`initial_subtract_operation` with
	``SUBTRACTED_DATA_ALL_NEW`` column, their dir-independent parmdbs, and the merged sky models in ``makesourcedb`` format.

Output
    For each band, all facet sources are added to the ``FACET_DATA_ALL`` column.

Pipeline Steps
    create_ms_map, create_parmdb_map, create_full_skymodels_map
        Make datamaps for input MS files, their dir-independent parmdbs, and
        the merged sky models.

    make_facet_skymodels_all
        Select all model components belonging to all facet sources and write
        these sky models in ``makesourcedb`` format.

    add_all_facet_sources
        Run BBS to add sources to the data using the sky models above.


Add sources from ``MODEL_DATA``
-------------------------------

.. note::

    This step is done only for the reimaging pipeline.

Input
	MS files from the :ref:`initial_subtract_operation` with ``SUBTRACTED_DATA_ALL_NEW`` column, their dir-dependent parmdbs, and the dataset with the self-calibration facet model phase shifted back to the field center.

Output
    For each band, all facet sources are added to the ``FACET_DATA_ALL`` column.

Pipeline Steps
    copy_column
        Copy the ``MODEL_DATA`` column from self calibration to the input datasets.

    expand_parmdb_map, add_all_facet_sources
        Run BBS to add sources to the data using the sky models above and the dir-dependent parmdbs.


Split and phase shift data
--------------------------

Input
	MS files from the previous step with ``FACET_DATA_ALL`` columns.

Output
    For each band, a phase-shifted dataset with all sources in the ``DATA`` column.

Pipeline Steps
    shift_all
        Run DPPP to split and phase shift ``FACET_DATA_ALL`` to the RA and Dec of the facet.


Imaging preparation
-------------------

Input
	MS files with phase-shifted facet sources in the ``DATA`` column and the dir-dependent parmdbs from self calibration.

Output
    Datasets ready for imaging.

Pipeline Steps
    expand_merged_parmdb_map, apply_dir_dep
        Apply dir-dependent solutions from self calibration to the phase-shifted ``DATA`` column to make a ``CORRECTED_DATA`` column for imaging.

    average, create_compressed_mapfile, concat_averaged
        Average ``CORRECTED_DATA`` column in time and frequency and concatenate in frequency in preparation for imaging.


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
        WSClean imaging run. Imaging is done with a cell size of 1.5". Wide-band imaging is done if more than 5 bands are used.


