Introduction
============

This document describes the facet calibration software that implements the facet calibration scheme described in van Weeren et al. (in prep.) and Williams et al. (in prep.).  The software are divided into a number of pipelines, each of which is described in the following sections.

In addition to these pipelines, the Factor software package has been developed to handle the setup and running of pipelines to perform facet calibration in a user-friendly way (note that currently the pre-facet calibration pipeline is not part of Factor). For instructions on using Factor, please see the Factor chapter in the LOFAR Imaging Cookbook.

The pre-facet calibration pipeline parset is available at https://github.com/tammojan/facet-calibration. The other pipeline parsets referenced here are available from the GitHub repository at https://github.com/revoltek/factor under the ``factor/pipeline/parsets`` directory.

.. note::
    If desired, Factor can be run with the ``-d`` option to generate pipeline parsets which can then be run by hand.

The overall structure of the pipeline after pre-facet calibration is shown in the Figure :num:`factor-flowchart` below. The pipeline is divided into a number of operations (subpipelines), the division of which is largely determined by whether or not multiple operations may be run in parallel. In this flowchart, each operation is outlined with a black box and is named with the corresponding pipeline parset name.

.. _factor-flowchart:

.. figure:: factor_flow.pdf
   :figwidth: 90 %
   :align: center

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
