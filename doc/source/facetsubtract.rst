.. _subtract_facet_sources:

Subtract Facet Sources Operation
================================

This operation subtracts all facet sources from the data. There are two possible pipeline parsets for this operation:

``facetsub_pipeline.parset``
    Full subtract operation. Facet sources are added to the empty data and then subtracted with the improved model and calibration from self calibration.

``facetsub_single_pipeline.parset``
    For use when the facet under consideration has been run singly (not in parallel) or is the first of a group of facets run in parallel. In this case, the full-res subtracted column from self calibration can be subtracted from the original data after phase shifting it back to the field center. For the other facets in a group, the facet sources have to be added and subtracted to pick up the improved subtracted data from the other directions in the group.

.. note::

    There should be one pipeline per facet, and each pipeline must be run in series as they each write data to the input datasets.

    This operation is separated from the self-calibration operation (:ref:`facet_selfcal`) as these pipelines
    must be run in series, whereas the self-calibration pipelines can be run in parallel.

Add sources
-----------

Input
	MS files from the :ref:`initial_subtract_operation` with
	``SUBTRACTED_DATA_ALL`` (or ``SUBTRACTED_DATA_ALL_NEW`` if at least one facet has gone through self calibration previously) columns, their dir-dependent self-calibration parmdbs, and the phase-shifted datasets with the FT-ed ``MODEL_DATA`` column.

Output
    For each band, all facet sources are added to the ``SUBTRACTED_DATA_ALL`` (or ``SUBTRACTED_DATA_ALL_NEW``) columns to create the ``FACET_DATA_ALL`` column.

Pipeline Steps
    update_shifted_all_hosts, update_input_bands_hosts
        Update hosts for various datamaps to reflect that the pipelines are now running in series.

    add_all_facet_sources
        Add all facet sources to the ``SUBTRACTED_DATA_ALL`` (or ``SUBTRACTED_DATA_ALL_NEW``) columns to make a ``FACET_DATA_ALL`` column.

Test data
    TODO


Subtract sources
----------------

Input
	MS files from the previous step with ``FACET_DATA_ALL`` column.

Output
    For each band, all facet sources are subtracted from the ``FACET_DATA_ALL`` column to produce a new ``SUBTRACTED_DATA_ALL_NEW`` column. Additionally, a dataset is made with the facet model phase shifted back to the field center (for use later if re-imaging of the facet is desired).

Pipeline Steps
    shift_to_field
        Phase shift the ``MODEL_DATA`` column from self calibration back to the field center.

    copy_column
        Copy the phase-shifted model column to the ``MODEL_DATA`` column of the input datasets.

    expand_parmdb_map
        Match the number of dir-dependent parmdb entries in the datamap to that in the input datasets datamap.

    subtract
        Subtract the ``MODEL_DATA`` column from the ``SUBTRACTED_DATA_ALL`` (or ``SUBTRACTED_DATA_ALL_NEW``) column.

    copy_shifted_map
        Copy datamap for shifted model dataset to convenient location.

Test data
    TODO


Partial subtract operation
--------------------------

.. note::

    This step is done only for the ``facetsub_single_pipeline.parset`` pipeline and replaces all of the above steps.

Input
	Phase-shifted datasets from the :ref:`facet_selfcal` with ``SUBTRACTED_DATA`` and ``MODEL_DATA`` columns.

Output
    For each band, a new ``SUBTRACTED_DATA_ALL_NEW`` column is made. Additionally, a dataset is made with the facet model phase shifted back to the field center (for use later if re-imaging of the facet is desired).

Pipeline Steps
    update_shifted_all_hosts, update_input_bands_hosts
        Update hosts for various datamaps to reflect that the pipelines are now running in series.

    shift_subdata_to_field
        Phase shift the ``SUBTRACTED_DATA`` column from self calibration back to the field center

    copy_column
        Copy the ``SUBTRACTED_DATA`` column to the ````SUBTRACTED_DATA_ALL_NEW```` column of the input datasets

    shift_model_to_field
        Phase shift the ``MODEL_DATA`` column from self calibration back to the field center

    copy_shifted_model_map
        Copy datamap for shifted model dataset to convenient location.

Test data
    TODO



