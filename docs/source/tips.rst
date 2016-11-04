.. _tips:

Tips for running Factor
=======================

Below are a number of tips for setting up and running Factor in the optimal way.

Match the number of chunks to the number of cores
-------------------------------------------------
To maximize CPU usage, the number of chunks should be evenly divisible by
the number of available cores. Adjust the :term:`chunk_size_sec` parameter to obtain the
desired number of chunks.

Use the ram drive
-----------------
Set :term:`dir_local_selfcal` = ``/dev/shm`` to use the ram drive, which will keep data used for self calibration in memory. Note that you must have sufficient space on ``/dev/shm`` to hold the data, which typically use ~ 5-10 GB per direction for a full-bandwidth run.

Use compression
---------------
Set :term:`use_compression` = ``True`` to reduce file sizes and IO with Dysco compression. Enabling this setting is particularly helpful for systems with slow IO.

Use baseline-dependent averaging in WSClean
-------------------------------------------
Set :term:`wsclean_bl_averaging` = ``True`` to speed up imaging by factors of 5-10.

Self calibrate in parallel
--------------------------

Self calibrate multiple directions at the same time by setting the groupings and ndir_per_node parameters appropriately for your system. Typically, ~ 8-16 cores per direction works well. Depending on the distribution of calibrator fluxes in your field, you may also be able to set the :term:`groupings` parameter to define groups that are several times larger than :term:`ndir_per_node` (or than the number of nodes if you are using a cluster). For example, for a field with only one very bright source (S > 5 Jy) and :term:`ndir_per_node` = 4, one can set :term:`groupings` = ``1:1,12:0``. This setup would result in the brightest source being processed and subtracted first and then the next 12 sources being subtracted only once all 12 have been processed.

Image only the target
---------------------

If you are interested only in a single target in the field, the following strategy is optimal:

    * Reorder the directions in the directions file so that the facets you wish to process are at the top. Generally, these facets should be those that contain the brightest sources and those that neighbor the target facet (or include the target). Note, however, that you do not need to put the target facet last since it will be automatically imaged after all directions have been processed.
    * Set the maximum number of selfcal cycles with :term:`max_selfcal_loops` to a low number (2 or 3) to speed up self calibration for non target sources. Spending a lot of time self calibrating non-target sources will generally not significantly improve the noise in the target facet (unless the calibrator is very bright and nearby). The number of selfcal cycles used for the target can be controlled with the :term:`target_max_selfcal_loops` option, which defaults to 10.
    * Set :term:`image_target_only` = ``True`` to skip imaging of the non-target facets.



