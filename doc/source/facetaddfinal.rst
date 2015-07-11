.. _add_final_facet_sources:

Add Final Facet Sources Operation
=================================

This operation adds the facet sources to the final empty data
in preparation for imaging (or re-imaging). There are two possible pipeline parsets for this operation:

``facetaddfinal_cc_skymodel_pipeline.parset``
    For use when the facet under consideration has not successfully gone through self calibration. In this case, the self-calibration solutions from the nearset facet should be used instead.

``facetaddfinal_model_image_pipeline.parset``
    For use when the facet under consideration has successfully gone through self calibration and reimaging is desired.

.. note::

    There should be one pipeline per facet, and each pipeline must be run in series as they each write data to the input datasets.

    This operation is separated from the final imaging operation (:ref:`facet_image`) as these pipelines
    must be run in series, whereas the imaging pipelines can be run in parallel.


Add sources from CC model
-------------------------

.. note::

    This step is done only for the ``facetaddfinal_cc_skymodel_pipeline.parset`` pipeline.

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

Test data
    With ``Test_data/RX42_SB070-079.2ch10s.ms``, this step produces the sky model ``NEP_SB070-079.2ch10s.wsclean_low2-model.make_facet_skymodels_all`` in ``Test_run/results/facetadd/facet_patch_***/``, which in turn is used to make the ``FACET_DATA_ALL`` column in this MS file.


Add sources from ``MODEL_DATA``
-------------------------------

.. note::

    This step is done only for the ``facetaddfinal_model_image_pipeline.parset`` pipeline.

Input
	MS files from the :ref:`initial_subtract_operation` with ``SUBTRACTED_DATA_ALL_NEW`` column, their dir-dependent parmdbs, and the dataset with the self-calibration facet model phase shifted back to the field center.

Output
    For each band, all facet sources are added to the ``FACET_DATA_ALL`` column.

Pipeline Steps
    copy_column
        Copy the ``MODEL_DATA`` column from self calibration to the input datasets.

    expand_parmdb_map, add_all_facet_sources
        Run BBS to add sources to the data using the sky models above and the dir-dependent parmdbs.

Test data
    With ``Test_data/RX42_SB070-079.2ch10s.ms`` and the phase-shifted facet all-source MS files (e.g., ``NEP_SB070-079.2ch10s.shift_all``) in ``Test_run/results/facetselfcal/facet_patch_***/``, this step produces the ``FACET_DATA_ALL`` column in the ``Test_data/RX42_SB070-079.2ch10s.ms`` file.


Split and phase shift data
--------------------------

Input
	MS files from the previous step with ``FACET_DATA_ALL`` columns.

Output
    For each band, a phase-shifted dataset with all no sources in the ``DATA`` column.

Pipeline Steps
    shift_all
        Run DPPP to split and phase shift ``FACET_DATA_ALL`` to the RA and Dec of the facet.

    copy_all_map
        Copy datamap for shifted datasets to convenient location.

Test data
    With ``Test_data/RX42_SB070-079.2ch10s.ms``, this step produces the MS files ``Test_run/results/facetaddfinal/facet_patch_***/NEP_SB070-079.2ch10s.shift_all``.

