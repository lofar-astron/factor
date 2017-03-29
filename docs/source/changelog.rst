.. _changelog:

Changelog
=========

Changes since version 1.1:

    * An unarchiving tool (``unarchivefactor``, see :ref:`unarchivefactor`) has been added that can unarchive an archive made with ``archivefactor``
    * An option (:term:`update_selfcal_clean_regions`) has been added that controls whether a user-supplied clean mask is updated during selfcal
    * Intersections in the user-supplied clean mask are now detected and an error raised
    * An archiving tool (``archivefactor``, see :ref:`archivefactor`) has been added that can archive the subtracted datasets, the sky models, the instrument tables, the selfcal plots, and the calibrated data for one or more directions
    * Polynomial sky models generated directly by WSClean during imaging are used for prediction, resulting in improved and faster subtraction of extended sources


Version 1.1
-------------

Changes since version 1.0:

    * The subtraction of sources after self calibration has been improved
    * The baseline-dependent averaging used during imaging has been reduced, as the previous averaging caused significant smearing away from the facet center
    * Many bug fixes


Version 1.0
-----------

First release. Changes since version pre1.0:

    * The BBS calibration software has been replaced by GainCal, which is both much faster and more stable
    * The CASA imager has been replaced by WSClean, which is generally much faster and integrates better with the LOFAR pipeline framework that Factor uses
    * A new option has been added to keep the primary data files used in self calibration in memory, speeding up self calibration dramatically on some systems (set :term:`dir_local_selfcal` = ``/dev/shm`` under the ``[cluster]`` section of the parset)
    * The Dysco storage manager can now be used to compress visibilities and weights, reducing file sizes (and IO) by a factor of ~ 2.5
    * The ``reimage_selfcaled`` option has been removed, as facets are always reimaged now. A :term:`image_target_only` option has been added if you wish to reimage only the target
    * The ``skip_facet_imaging`` option has been removed, as facets are always imaged now to improve the subtraction (using only a faction of the full bandwidth)

