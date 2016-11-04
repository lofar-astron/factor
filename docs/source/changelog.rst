.. _changelog:

Changelog
=========

Version 1.0
-----------

First release. Changes since version pre1.0:

    * The BBS calibration software has been replaced by GainCal, which is both much faster and more stable
    * The CASA imager has been replaced by WSClean, which is generally much faster and integrates better with the LOFAR pipeline framework that Factor uses
    * A new option has been added to keep the primary data files used in self calibration in memory, speeding up self calibration dramatically on some systems (set :term:`dir_local_selfcal` ``= /dev/shm`` under the ``[cluster]` section of the parset)
    * The Dysco storage manager can now be used to compress visibilities and weights, reducing file sizes (and IO) by a factor of ~ 2.5
    * The ``reimage_selfcaled`` option has been removed, as facets are always reimaged now. A :term:`image_target_only` option has been added if you wish to reimage only the target
    * The ``skip_facet_imaging`` option has been removed, as facets are always imaged now to improve the subtraction (using up to 6 bands)

