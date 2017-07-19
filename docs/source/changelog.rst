.. _changelog:

Changelog
=========

Changes since version 1.3:

    * Options have been added to allow the scales used by WSClean during clean to be specified (selfcal_multiscale_scales_pixel and facet_multiscale_scales_pixel under [imaging])
    * The atrous_do column in the directions file has been renamed to mscale_selfcal_do, as this name better describes its purpose
    * Improved handling of flagged solutions during smoothing of the amplitude solutions


Version 1.3
-----------

Changes since version 1.2:

    * Updated to use WSClean v2.4 and LSMTool v1.2.0
    * Preaveraging is now done in frequency as well as in time. This preaveraging generally improves the S/N of the CS slow-gain solutions for fainter sources by a factor of ~ 2
    * An option (:term:`min_fraction_per_band`) has been added that sets the minimum allowed unflagged fraction per band
    * WSClean's automasking feature is now used during imaging. The old image-mask-image sequence is no longer used during self calibration, but can still be used during the final, full-bandwidth facet imaging if `automask_facet_image = False` under the `[imaging]` section of the parset
    * Handling of pipeline failure/interruption has been improved
    * The combination of flagging ranges specified by the `flag_abstime`, `flag_baseline`, and `flag_freqrange` options can now be set with the `flag_expr` option


Version 1.2
-----------

Changes since version 1.1:

    * The combination of flagging ranges specified by the :term:`flag_abstime`, :term:`flag_baseline`, and :term:`flag_freqrange` options can now be set with the :term:`flag_expr` option
    * An unarchiving tool (``unarchivefactor``, see :ref:`unarchivefactor`) has been added that can unarchive an archive made with ``archivefactor``
    * An option (`update_selfcal_clean_regions``) has been added that controls whether user-supplied clean masks are updated during selfcal
    * Intersections in user-supplied clean masks are now detected and an error raised
    * An archiving tool (``archivefactor``, see :ref:`archivefactor`) has been added that can archive the subtracted datasets, the sky models, the instrument tables, the selfcal plots, and the calibrated data for one or more directions
    * Polynomial sky models generated directly by WSClean during imaging are used for prediction, resulting in improved and faster subtraction of extended sources


Version 1.1
-----------

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

