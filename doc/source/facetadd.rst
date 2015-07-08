.. _add_facet_sources:

Add Facet Sources Operation
===========================

This operation adds both the facet calibrator and all facet sources to the data
in preparation for self calibration. The pipeline parset for this operation is ``facetadd_pipeline.parset``.

.. note::

    There should be one pipeline per facet, and each pipeline must be run in series as they each write data to the input datasets.

    This operation is separated from the self-calibration operation (:ref:`facet_selfcal`) as these pipelines
    must be run in series, whereas the self-calibration pipelines can be run in parallel.

Add sources
-----------

Input
	MS files from the :ref:`initial_subtract_operation` with
	``SUBTRACTED_DATA_ALL`` (or ``SUBTRACTED_DATA_ALL_NEW`` if at least one facet has gone through self calibration previously) and ``CORRECTED_DATA`` columns, their dir-independent parmdbs, and the merged sky models in ``makesourcedb`` format.

Output
    For each band, the facet calibrator and all facet sources are added to the ``FACET_DATA_CAL`` and ``FACET_DATA_ALL`` columns, respectively.

Pipeline Steps
    create_ms_map, create_parmdb_map, create_full_skymodels_map
        Make datamaps for input MS files, their dir-independent parmdbs, and
        the merged sky models.

    make_facet_skymodels_all, make_facet_skymodels_cal
        Select all model components belonging to the facet calibrator and to all facet sources and writes
        these sky models in ``makesourcedb`` format.

    add_all_facet_sources, add_cal_facet_sources
        Run BBS to add sources to the data using the sky models above.

Test data
    TODO


Split and phase shift data
--------------------------

Input
	MS files from the previous step with ``FACET_DATA_CAL`` and ``FACET_DATA_ALL`` columns.

Output
    For each band, three phase-shifted datasets with calibrator, all, and no sources in the ``DATA`` column.

Pipeline Steps
    shift_all, shift_cal, shift_empty
        Run DPPP to split and phase shift ``FACET_DATA_CAL``, ``FACET_DATA_ALL``, and
        ``SUBTRACTED_DATA_ALL`` columns (or ``SUBTRACTED_DATA_ALL_NEW`` column if at least one facet has gone through self calibration previously) to the RA and Dec of the facet.

    copy_all_map, copy_cal_map, copy_empty_map
        Copy datamaps for shifted datasets to convenient location.

Test data
    TODO


