Introduction
============

The facet calibration software is divided into a number of pipelines, each of which is described in the following sections. The pre-facet calibration pipeline parset is available at https://github.com/tammojan/facet-calibration. The other pipeline parsets referenced here are available from the GitHub repository at https://github.com/revoltek/factor under the ``factor/pipeline/parsets`` directory.

In addition to these pipelines, the Factor software package has been developed to handle the setting up and running of pipelines to perform facet calibration in a user-friendly way (note that currently the pre-facet calibration pipeline is not part of Factor). For instructions on using Factor, please see the sections below and the Factor chapter in the LOFAR Imaging Cookbook.

The overall structure of the pipeline after pre-facet calibration is shown in the Figure :num:`factor-flowchart` below. The pipeline is divided into a number of operations (subpipelines), the division of which is largely determined by whether or not multiple operations may be run in parallel. In this flowchart, each operation is outlined with a black box and is named with the corresponding pipeline parset name.

.. _factor-flowchart:

.. figure:: factor_flow.pdf
   :figwidth: 90 %
   :align: center

   Factor flowchart


The Factor Parset
-----------------

Before Factor can be run, a parset describing the reduction should be made. An example parset for a basic reduction on a single machine is shown below::

        [global]
        dir_working = /path/to/factor/working/dir
        dir_ms = /path/to/prefacet/working/dir

        [directions]
        flux_min_jy = 0.1
        size_max_arcmin = 3.0
        separation_max_arcmin = 7.5
        max_num = 30
        groupings = 2:0

        [cluster]
        ndir_per_node = 2

In the global section of this parset, the ``dir_ms`` directory should be the directory containing the output of the pre-facet preparation pipeline. Factor will automatically identify the final MS files from the pre-facet preparation (those with the extension ``.dpppconcat``). In the directions section, various parameters governing how DDE calibrators are identified are set. Here, calibrators must meet the following criteria: a minimum total flux density of 0.1 Jy and a maximum size of 3 arcmin, both in the highest-frequency band. Calibrators within 7.5 arcmin of each other will be combined together into a single facet. The maximum number of facets is set to 30, and finally the facets are grouped into groups of 2 for self calibration. The facets in a group are self calibrated simultaneously, speeding up processing if enough cores and memory are available. However, the source-subtracted data for the facets in a group are identical, and hence the artifacts from one facet in the group may affect the calibration of the other facets. To minimize this possibility, reordering is done to ensure that the facets in a group lie far from one another.

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


Running Factor
--------------

Factor can now be run with::

    $ runfactor factor.parset

where ``factor.parset`` is the parset discussed above.

.. note::
    If desired, Factor can be run with the ``-d`` option to generate pipeline parsets which can then be run by hand.

Factor will begin by checking the input MS files and the direction-independent instrument tables. If the instrument tables contain real/imaginary values, they are converted the phase/amplitude. Next, Factor will check whether the initial subtract operation (described in detail in :ref:`initial_subtract_operation`) is required and will run it if needed.

Next, Factor will check for a file describing the DDE calibrators. If not found, Factor will select the DDE calibrators automatically and generate the facet regions. At this point, if ``interactive = True`` is set in the parset, Factor will pause to allow a check on the facets.

Factor will then begin self calibration and imaging (described in detail in :ref:`facet_selfcal`) of the first facet or group of facets. For the facet or facets that pass self-calibration verification, Factor subtracts the improved model using the direction-dependent instrument tables (described in detail in :ref:`subtract_facet_sources`). Processing then proceeds to the self calibration and imaging of the next facet or facet group, and the step are looped until all facets have been processed.

After self calibration is finished, imaging is done for any facets that did not successfully go through self calibration. These facets receive the direction-dependent instrument tables of the nearest facet for which self calibration succeeded. Lastly, all facet images are mosaicked together and the primary beam attenuation is removed to produce the final image of the entire field.
