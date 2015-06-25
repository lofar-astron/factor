.. _facet_selfcal:

Facet Self-Calibration Operation
================================

This operation performs self calibration on the facet calibrator.

.. note::

    There should be one pipeline per facet, and the pipelines may be run in parallel.

    This operation is separated from the add (:ref:`add_facet_sources`) and subtract (:ref:`subtract_facet_sources`) operations as those pipelines
    must be run in series.


.. _selfcal_cycle:

Self-calibration cycle
----------------------

Input
	MS files from the :ref:`initial_subtract_operation` with
	``SUBTRACTED_DATA_ALL`` (or ``SUBTRACTED_DATA_ALL_NEW`` if at least one facet has gone through self calibration previously) and ``CORRECTED_DATA`` columns, their dir-independent parmdbs, and the merged sky models in ``makesourcedb`` format.

Output
    For each band, the facet calibrator and all facet sources are added to the ``FACET_DATA_CAL`` and ``FACET_DATA_ALL`` columns, respectively.

Pipeline Steps
    create_ms_map, create_parmdb_map, create_full_skymodels_map
        Make datamaps for input MS files, their dir-independent parmdbs, and
        the merged sky models

    make_facet_skymodels_all, make_facet_skymodels_cal
        Selects all model components belonging to the facet calibrator and to all facet sources and writes
        these sky models in ``makesourcedb`` format.

    add_all_facet_sources, add_cal_facet_sources
        Run BBS to add sources to the data using the sky models above.

Test data
    TODO

