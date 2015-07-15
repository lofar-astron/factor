
Factor: Facet Calibration Pipeline
==================================

This document describes the Factor facet calibration pipeline that implements the facet calibration scheme described in van Weeren et al. (2015). All of pipeline parsets referenced here are available from the GitHub repository at https://github.com/revoltek/factor under the ``factor/pipeline/parsets`` directory. The pipeline is divided into a number of operations (subpipelines), the division of which is largely determined by whether or not multiple operations may be run in parallel.

The overall structure of the pipeline is shown in the `Factor flowchart`_ below. In this flowchart, each operation is outlined with a black box and is named with the corresponding pipeline parset name.

.. _`Factor flowchart`:

.. figure:: factor_flow.pdf
   :scale: 40 %
   :figwidth: 75 %
   :align: center
   :alt: example image

   Factor flowchart

.. note::

    The test data (and the resulting outputs) referenced in this document are available at ``lof021:/data/scratch/rafferty/Factor_test``. The outputs can be reproduced with the following Factor parset::

        [global]
        dir_working = /data/scratch/rafferty/Factor_test/Test_run
        dir_ms = /data/scratch/rafferty/Factor_test/Test_data
        parmdb_name = instrument_ap_smoothed

        [directions]
        ndir_total = 3
        ndir_selfcal = 2
        flux_min_jy = 0.1
        size_max_arcmin = 3.0
        separation_max_arcmin = 7.5
        max_num = 100
        groupings = 2:0

        [cluster]
        ndir_per_node = 2

.. toctree::
   :maxdepth: 2
   :numbered:

   initsubtract.rst
   facetadd.rst
   facetselfcal.rst
   facetsubtract.rst
   facetaddfinal.rst
   facetimagefinal.rst
   mosaic.rst
   factor.rst

