
Factor: Facet Calibration Pipeline
==================================

This document describes each step in the Factor facet calibration pipeline. The pipeline is divided into a number of operations. Each operation is a group of steps. This grouping is determined by whether or not multiple pipelines may be run in parallel.

.. note::

    The test data referenced in this document is available at ``lof021:/data/scratch/rafferty/Factor_test``. The outputs can be reproduced with the following Factor parset::

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


Index and search
================

* :ref:`genindex`
* :ref:`search`

