.. _runfactor:

Running the Factor Software
===========================

Once the pre-facet preparation (see Section :ref:`pre_facet`) is complete, the Factor software package can be used to perform the facet calibration and imaging. The following sections describe the steps needed before Factor can be run.


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

In the global section of this parset, the ``dir_ms`` directory should be the directory containing the output of the pre-facet preparation pipeline. Factor will automatically identify the final MS files from the pre-facet preparation (those with the extension ``.dpppconcat``). In the directions section, various parameters governing how DDE calibrators are identified are set. Here, calibrators must meet the following criteria: a minimum total flux of 0.1 Jy in the highest-frequency band and a maximum size of 3 arcmin. Calibrators withing 7.5 arcmin of each other will be combined together into a single facet. The maximum number of facets is set to 30, and finally the facets are grouped into groups of 2 for self calibration. The facets in a group are self calibrated simultaneously, speeding up processing if enough cores and memory are available but with


Running Factor
--------------

Factor can now be run with::

    $ runfactor factor.parset
