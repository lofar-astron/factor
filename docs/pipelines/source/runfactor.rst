.. _runfactor:

Running the Factor Software
===========================

Once the pre-facet preparation (see Section :ref:`pre_facet`) is complete, the Factor software package can be used to perform the facet calibration and imaging. The following sections describe the steps needed before Factor can be run.


The Factor Parset
-----------------

Before Factor can be run, a parset describing the reduction should be made. An example parset is shown below::

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

Factor will automatically identify the final MS files from the pre-facet preparation (those with the extension ``.dpppconcat``).


Running Factor
--------------

Factor can now be run with::

    $ runfactor factor.parset
