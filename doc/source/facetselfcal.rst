.. _facet_selfcal:

Facet Self-Calibration Operation
================================

This operation performs self calibration on the facet calibrator.

.. note::

    There should be one pipeline per facet, and the pipelines may be run in parallel.

    This operation is separated from the add (:ref:`add_facet_sources`) and subtract (:ref:`subtract_facet_sources`) operations as those pipelines
    must be run in series.


Data preparation
----------------

Input
	MS files from the :ref:`facet_add` operation with phase-shifted facet calibrator in the
	``DATA`` column and the dir-independent parmdbs.

Output
    Datasets ready for self calibration.

Pipeline Steps
    apply_dir_indep
        Apply dir-independent solutions to the phase-shifted ``DATA`` column to make a ``CORRECTED_DATA`` column for imaging.

    average_data, average_corr
        Average ``DATA`` and ``CORRECTED_DATA`` columns to 1 channel per band.

    create_compressed_mapfile_data, create_compressed_mapfile_corr
        Create datamaps suitable for DPPP concatenation.

    concat_data, concat_corr
        Run DPPP to concatenate all bands together

    copy_column
        Copy the ``CORRECTED_DATA`` column so that a single MS file has all needed columns.

Test data
    TODO


.. _selfcal_cycle:

Self-calibration cycle
----------------------
The general self-calibration cycle is described here. Modification to this cycle
are described in later steps.

Input
	Concatenated, averaged MS file with ``DATA`` and ``CORRECTED_DATA`` columns

Output
    Improved ``MODEL_DATA`` column and dir-dependent solutions.

Pipeline Steps
    averageX
        Average in time in preparation for imaging

    casa_imageX1, adjust_casa_mapfileX, maskX, casa_imageX2
        CASA imaging run. Imaging is done with a cell size of 1.5". Wide-band imaging is done if more than 5 bands are used. Multi-scale clean is always used.

    create_modelX_map, casa_ftX
        CASA FT run

    solveX
        Solve for phases or amplitudes


Self-calibration cycle 0
------------------------
The self-calibration cycle (see :ref:`selfcal_cycle`) is performed with phase-only calibration. The resulting image should be similar to the dir-independent images obtained in the :ref:`initial_subtract_operation` (although of higher resolution). An example image is shown in the `Cycle 0 example image`_.

.. _`Cycle 0 example image`:

.. figure:: cycle_0_image.png
   :scale: 40 %
   :figwidth: 75 %
   :align: center
   :alt: example image

   Cycle 0 example image


Self-calibration cycle 1
------------------------
The self-calibration cycle (see :ref:`selfcal_cycle`) is performed with phase-only calibration. The resulting image should show marked improvement over the cycle-0 image. An example image is shown in the `Cycle 1 example image`_.

.. _`Cycle 1 example image`:

.. figure:: cycle_1_image.png
   :scale: 40 %
   :figwidth: 75 %
   :align: center
   :alt: example image

   Cycle 1 example image


Self-calibration cycle 2
------------------------
The self-calibration cycle (see :ref:`selfcal_cycle`) is performed with phase-only calibration. The resulting image may or may not show improvement over the cycle-1 image. An example image is shown in the `Cycle 2 example image`_.

.. _`Cycle 2 example image`:

.. figure:: cycle_2_image.png
   :scale: 40 %
   :figwidth: 75 %
   :align: center
   :alt: example image

   Cycle 2 example image


Self-calibration cycle 3
------------------------
The self-calibration cycle (see :ref:`selfcal_cycle`) is performed with phase and amplitude calibration (fast phase, slow amplitude). The resulting image should show marked improvement over the cycle-2 image. An example image is shown in the `Cycle 3 example image`_.

.. note::

    Negative features in the image are due to poorly subtracted sources from the :ref:`initial_subtract_operation`.

.. _`Cycle 3 example image`:

.. figure:: cycle_3_image.png
   :scale: 40 %
   :figwidth: 75 %
   :align: center
   :alt: example image

   Cycle 3 example image


Smooth amplitudes 1
-------------------
The slow amplitude solutions from cycle 3 are smoothed to remove outliers.


Self-calibration cycle 4
------------------------
The self-calibration cycle (see :ref:`selfcal_cycle`) is performed with phase and amplitude calibration (fast phase, slow amplitude). The resulting image may or may not show improvement over the cycle-3 image. An example image is shown in the `Cycle 4 example image`_.

.. _`Cycle 4 example image`:

.. figure:: cycle_4_image.png
   :scale: 40 %
   :figwidth: 75 %
   :align: center
   :alt: example image

   Cycle 4 example image


Smooth amplitudes 2
-------------------
The slow amplitude solutions from cycle 4 are smoothed to remove outliers.


Merge self-calibration parmdbs
------------------------------

Input
	Fast phase and slow amplitude solution parmdbs.

Output
    Merged parmdb with both fast phase and slow amplitude solutions.

Pipeline Steps
    merge_selfcal_parmdbs
        Apply dir-independent solutions to the phase-shifted ``DATA`` column to make a ``CORRECTED_DATA`` column for imaging. An example of the solutions for RS106 is shown in `Merged parmdb fast solutions plot`_ and `Merged parmdb slow solutions plot`_.

    .. _`Merged parmdb fast solutions plot`:

    .. figure:: merged_parmdb_fast_plot.png
       :scale: 40 %
       :figwidth: 75 %
       :align: center
       :alt: example solutions

       Merged parmdb fast solutions plot

    .. _`Merged parmdb slow solutions plot`:

    .. figure:: merged_parmdb_slow_plot.png
       :scale: 40 %
       :figwidth: 75 %
       :align: center
       :alt: example solutions

       Merged parmdb slow solutions plot


Test data
    TODO
