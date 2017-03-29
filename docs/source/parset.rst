.. _factor_parset:

The Factor parset
=================

Before Factor can be run, a parset describing the reduction must be made. The
parset is a simple text file defining the parameters of a run in a number of
sections. For example, a typical parset for a basic reduction on a single
machine could look like the following (see :ref:`tips` for tips on setting up an
optimal parset):

.. code-block:: none

        [global]
        dir_working = /path/to/factor/working/dir
        dir_ms = /path/to/prefacet/working/dir

        [directions]
        flux_min_jy = 0.3
        size_max_arcmin = 1.0
        separation_max_arcmin = 9.0
        max_num = 40

The available options are described below under their respective sections.


.. _parset_global_options:

``[global]``
------------

.. glossary::

    dir_working
        Full path to working dir where Factor will run (required). All output
        will be placed in this directory. E.g., ``dir_working = /data/wdir``.

    dir_ms
        Full path to directory containing input bands. It will be scanned for
        all ``.MS`` and ``.ms`` files (required). Note that FACTOR works on a
        copy of these files and does not modify the originals in any way. E.g.,
        ``dir_ms = /data/bands``.

    parmdb_name
        Parmdb name for dir-indep. selfcal solutions (stored inside the input
        band measurement sets, so path should be relative to those; default =
        ``instrument_directionindependent``).

    skymodel_extension
        Extension that when concatenated with the "extension-stripped" MS path gives
        a path that is checked if it contains a skymodel. The default finds the skymodel
        files from the standard prefactor ``Initial-Subtract.parset``
        (default = ``.wsclean_low2-model.merge`` ; note the leading ".").

    chunk_size_sec
        Size of time chunks in seconds (default = 2400; minimum allowed value is
        1200). Ideally, the number of chunks should be evenly divisible by the
        total number of CPUs available to each direction (controlled by the
        options under :ref:`parset_cluster_options`). To prevent Factor from
        chunking the data, set this value to be larger than the length of the
        longest dataset (in this case, Factor will not make copies of the files
        but will make symbolic links to them instead, so please make backup
        copies yourself).

    use_compression
        Use Dysco compression for chunked files (default = False). Enabling this
        option will result in less storage usage and signifcanctly faster
        processing on systems with slow IO. To use this option, you must have the
        Dysco library in your ``LD_LIBRARY_PATH``. Note: if enabled, Factor will not
        make symbolic links to the input data, even if they are shorter than
        :term:`chunk_size_sec`, but will copy them instead.

    interactive
        Use interactive mode (default = ``False``). If ``True``, Factor will ask for confirmation of
        internally derived DDE calibrators and facets.

    keep_avg_facet_data
        Keep averaged calibrated data for each facet to allow re-imaging by hand
        (default = ``True``). If a target is specified (see below), the averaged
        data for the target is always kept, regardless of this setting. If the
        averaged data are kept, reimaging will be dramatically faster if
        multiple images per facet are made (e.g., at different scales)

    keep_unavg_facet_data
        Keep unaveraged calibrated data for each facet (default = ``False``).

    flag_abstime
        Range of times to flag (default = no flagging). The syntax is that of
        the preflagger ``abstime`` parameter (see the DPPP documentation on the LOFAR wiki for
        details of the syntax). E.g., ``[12-Mar-2010/11:31:00.0..12-Mar-2010/11:50:00.0]``.
        Note that time, frequency (set with
        :term:`flag_freqrange`), and baseline (set with
        :term:`flag_baseline`) ranges are AND-ed to produce the final flags.

    flag_baseline
        Range of baselines to flag (default = no flagging). The syntax is that
        of the preflagger ``baseline`` parameter (see the DPPP documentation for
        details of the syntax). E.g., ``flag_baseline = [CS013HBA*]``. Note that
        baseline, frequency (set with
        :term:`flag_freqrange`), and time (set with
        :term:`flag_abstime`) ranges are AND-ed to produce the final flags.

    flag_freqrange
        Range of frequencies to flag (default = no flagging). The syntax is that
        of the preflagger ``freqrange`` parameter (see the DPPP documentation for
        details of the syntax). E.g., ``flag_freqrange = [125.2..126.4MHz]``. Note that
        frequency, baseline (set with
        :term:`flag_baseline`), and time (set with
        :term:`flag_abstime`) ranges are AND-ed to produce the final flags.


.. _parset_calibration_options:

``[calibration]``
-----------------

.. glossary::

    exit_on_selfcal_failure
        Exit if selfcal fails for any direction (default = ``True``). If ``False``, processing
        will continue and the failed direction will receive the selfcal solutions of
        the nearest successful direction.

    skip_selfcal_check
        Skip self calibration check (default = ``False``). If ``True``,
        processing continues as if the selfcal succeeded.

    max_selfcal_loops
        Maximum number of cycles of the last step of selfcal to perform (default =
        10). The last step is looped until the number of cycles reaches this value or
        until the improvement in dynamic range over the previous image is less than
        1.25%.

    target_max_selfcal_loops
        Maximum number of cycles of the last step of selfcal to perform for the target
        facet, if any (default = 10).

    preapply_first_cal_phases
        Preapply the direction-dependent phase solutions for the first calibrator to
        all subsequent ones (default = ``False``). If ``True``, residual clock errors are
        removed before calibration and a single TEC+CommonScalarPhase fit is used
        across the whole bandwidth.

    preaverage_flux_Jy
        Use baseline-dependent preaveraging to increase the signal-to-noise of the
        phase-only solve for sources below this flux density (default = 0.0; i.e.,
        disabled). When activated, averaging in time is done to exploit the time
        coherence in the TEC solutions.

    multires_selfcal
        Use multi-resolution selfcal that starts at 20 arcsec resolution and increases the
        resolution in stages to the full resolution (default = ``False``). This method may
        improve convergence, especially when the starting model is poor.

    TEC_block_MHz
        Size of frequency block in MHz over which a single TEC+CommonScalarPhase solution is fit
        (default = 10.0).

    peel_flux_Jy
        Peel the calibrator for sources above this flux density in Jy (default = 25.0).
        When activated, the calibrator is peeled using a supplied sky model and
        the facet is then imaged as normal. Note: for each source that should be
        peeled, a sky model must be specified in the directions file in the
        :term:`peel_skymodel` column or be one of those included in Factor; if not, the
        calibrator will go through self calibration as if it were a normal calibrator.

    solve_min_uv_lambda
        Minimum uv distance in lambda for calibration (default = 80.0).

    spline_smooth2D
        Smooth amplitudes with spline fit + 2-D median (default = ``True``). If
        ``False``, smoothing is done with a 1-D median.

    solve_all_correlations_flux_Jy
        Include XY and YX correlations during the slow gain solve for sources above
        this flux density (default = 1000.0; i.e., effectively off). Below this value,
        only the XX and YY correlations are included. Note that :term:`spline_smooth2D` must
        be ``True`` to solve for all correlations. If you want to use it, then an useful
        value would be, e.g., 5.0.


.. _parset_imaging_options:

``[imaging]``
-----------------

.. glossary::

    make_mosaic
        Make final mosaic (default = ``True``).

    image_target_only
        Image only the target facet (default = ``False``). If ``True`` and a target is
        specified in the :ref:`parset_directions_options` section, then only the facet containing the
        target source is imaged.

    wsclean_image_padding
        Padding factor for WSClean images (default = 1.6).

    max_peak_smearing
        Max desired peak flux density reduction at center of the facet edges due to
        bandwidth smearing (at the mean frequency) and time smearing (default = 0.15 =
        15% reduction in peak flux). Higher values result in shorter run times but
        more smearing away from the facet centers. This value only applies to the
        facet imaging (self calibration always uses a value of 0.15).

    wsclean_nchannels_factor
        Max factor used to set the number of WSClean channel images when wide-band
        clean is used (default = 4). The number of channel images is determined by
        dividing the number of bands by the nearest divisor to this factor. Smaller
        values produce better results but require longer run times. Wide-band clean is
        activated when there are more than 5 bands.

    fractional_bandwidth_selfcal_facet_image
        Fractional of bandwidth to use for facet imaging during selfcal (default =
        0.25). Facet imaging during selfcal is used to improve the subtraction of
        non-calibrator sources in the facet. More bandwidth will result in a better
        subtraction but also longer runtimes

    wsclean_bl_averaging
        Use baseline-dependent averaging in WSClean (default = ``True``). If enabled,
        this option can dramatically speed up imaging with WSClean.
        NOTE: this option requires WSClean v2.0 or higher.

    selfcal_cellsize_arcsec
        Self calibration pixel size in arcsec (default = 1.5).

    selfcal_robust
        Self calibration Briggs robust parameter (default = -0.5).

    selfcal_min_uv_lambda
        Self calibration minimum uv distance in lambda (default = 80).

    selfcal_clean_threshold
        Use a clean threshold during selfcal imaging (default = ``False``). If ``False``,
        clean will always stop at 1000 iterations. If ``True``, clean will stop when it
        reaches the 1 sigma noise level.

    selfcal_adaptive_threshold
        Use an adaptive masking threshold during selfcal imaging (default = ``False``). If
        ``True``, the masking threshold will be estimated using the negative peaks in the
        image, which can help selfcal convergence in the presence of strong artifacts.

    update_selfcal_clean_regions
        Update user-supplied clean regions (i.e., those specified in the
        directions file under the :term:`region_selfcal` column) with new
        regions found by the source finder during selfcal (default = ``True``) .
        Facet regions (specified in the term:`region_facet` column of the
        directions file) are always updated

.. note::

    The following four parameters can be specified as lists if more than one set
    of images is desired. In this case, they must all have the same number of
    entries.

    facet_cellsize_arcsec
        Facet image pixel size in arcsec (default = self calibration value). E.g.,
        ``facet_cellsize_arcsec = [1.5, 15.0]``.

    facet_robust
        Facet image Briggs robust parameter (default = self calibration value). E.g.,
        ``facet_robust = [-0.25, 0.0]``.

    facet_taper_arcsec
        Facet image uv taper in arcsec (default = self calibration value). E.g.,
        ``facet_taper_arcsec = [0.0, 45.0]``.

    facet_min_uv_lambda
        Facet image minimum uv distance in lambda (default = self calibration value). E.g.,
        ``facet_min_uv_lambda = [80.0, 160.0]``.


.. _parset_directions_options:

``[directions]``
-----------------

.. glossary::

    faceting_skymodel
        Full path to sky model (in makesourcedb format) to be used for calibrator
        selection and facet-boundary source avoidance (default is to use
        direction-independent sky model of the highest-frequency band). The sky
        model must be grouped into patches by source (in PyBDSM, this grouping can be
        done by setting ``bbs_patches = 'source'`` in the ``write_catalog`` task)

    max_radius_deg
        Radius from phase center within which to consider sources as potential
        calibrators (default = 2 * FWHM of primary beam of highest-frequency band).

    directions_file
        Full path to file containing calibrator directions. If not given, directions
        are selected internally using the flux density and size cuts that follow.

    flux_min_for_merging_Jy
        Minimum flux density in Jy of a source to be considered for merging with a
        nearby source to form a calibrator group (default = 0.1).

    separation_max_arcmin
        Maximum separation between sources in arcmin below which they are
        grouped into a calibrator group (no default).

    size_max_arcmin
        Maximum size of individual sources to be considered for grouping into a
        calibrator group (no default).

    flux_min_Jy
        Minimum total flux density of a source (or group) to be considered as a calibrator (no default).

    minimize_nonuniformity
        When identifying calibrators with the above selection criteria, search for the
        set of calibrators that minimizes non-uniformity (default = ``False``). Generally,
        enabling this option will result in facets that are more uniform in size

    ndir_max
        Number of internally derived directions can be limited to a maximum number
        of directions if desired (default = all).

    ndir_process
        Total number of directions to process (default = all). If this number is
        greater than :term:`ndir_selfcal`, then the remaining directions will not be selfcal-
        ed but will instead be imaged with the selfcal solutions from the nearest
        direction for which selfcal succeeded (if a target is specified and
        :term:`target_has_own_facet` is ``True``, it will be imaged in this way after ndir_total
        number of directions are processed).

    ndir_selfcal
        Total number of directions to selfcal (default = all).

    faceting_radius_deg
        Radius within which facets will be used (default = 1.25 * FWHM / 2 of primary beam
        of highest-frequency band); outside of this radius, small patches are used
        that do not appear in the final mosaic.

    check_edges
        Check whether any sources from the initial subtract sky model fall on facet
        edges. If any are found, the facet regions are adjusted to avoid them (default
        is ``True``).

    groupings
        Grouping of directions into groups that are selfcal-ed in parallel, defined as
        grouping:n_total_per_grouping. For example, ``groupings = 1:5, 4:0`` means two
        groupings are used, with the first 5 directions put into groups of one (i.e.,
        each direction processed in series) and the rest of the directions divided
        into groups of 4 (i.e., 4 directions processed in parallel). Default is one at
        a time (i.e., ``groupings = 1:0``).

    allow_reordering
        If groups are used to process more than one direction in parallel, reordering
        of the directions in the groups can be done to maximize the flux-weighted
        separation between directions in each group (default = ``True``). This
        sorting attempts to minimize the effects that any artifacts from one
        direction might have on the other simultaneously processed directions.

    target_ra
        RA of the center of a circular region that encloses the target source
        (to ensure that it falls entirely within a single facet; no default). E.g.,
        ``target_ra = 14h41m01.884``.

    target_dec
        Dec of the center of a circular region that encloses the target source
        (to ensure that it falls entirely within a single facet; no default). E.g.,
        ``target_dec = +35d30m31.52``.

    target_radius_arcmin
        Radius in arcmin of a circular region that encloses the target source (to ensure
        that it falls entirely within a single facet; no default). Note that :term:`check_edges`
        must be True for the facet boundaries to be adjusted.

    target_has_own_facet
        The target can be placed in a facet of its own. In this case, it will
        not go through selfcal but will instead use the selfcal solutions of the
        nearest facet for which selfcal was done (default = ``False``).


.. _parset_cluster_options:

``[cluster]``
-----------------

.. glossary::

    clusterdesc_file
        Full path to cluster description file. Use ``clusterdesc_file = PBS`` to use the
        PBS / torque reserved nodes, clusterdesc_file = SLURM to use SLURM reserved
        ones, or use ``clusterdesc_file = JUROPA_slurm`` to use
        multiple nodes in a slurm reservation on JUROPA.
        If not given, the clusterdesc file for a single (i.e., local) node is used.

        .. note::

            On a cluster that uses PBS or SLRUM, Factor will automatically determine the nodes for which you have a
            reservation and use them. Note that you must ask for all the nodes you need
            in a single PBS or SLURM script, so that all nodes are available for the full Factor run. An
            example PBS script that uses 6 nodes (with 6 CPUs each) is shown below::

                #!/bin/bash
                #PBS -N Factor
                #PBS -l walltime=100:00:00
                #PBS -l nodes=6:ppn=6

                cd $PBS_O_WORKDIR
                source ~rafferty/init_factor
                runfactor factor.parset

    dir_local
        Full path to a local disk on the nodes for IO-intensive processing. The path
        must be the same for all nodes. Note: do not specify this parameter if you are
        running more than one direction simultaneously on a single machine, as it will cause conflicts between directions
        that are processed in parallel (no default).

    dir_local_selfcal
        Full path to ram drive (e.g., /dev/shm) to allow certain selfcal data to
        be cached in memory, speeding up selfcal on most systems considerably.

    ncpu
        Maximum number of CPUs per node to use (default = all). Note that this
        number will be divided among the directions to be run in parallel on
        each node (controlled by the :term:`ndir_per_node` option). Ideally, the
        number of time chunks (controlled by the :term:`chunk_size_sec` option)
        should be evenly divisible by the number of CPUs per direction.

    nthreads_io
        Maximum number of IO-intensive threads to run per node (default =
        sqrt(:term:`ncpu`)). Note that this number will be divided among the
        directions to be run in parallel on each node (controlled by the
        :term:`ndir_per_node` option). Ideally, the number of time chunks (controlled
        by the :term:`chunk_size_sec` option) should be evenly divisible by the
        number of IO-intensive threads per direction.

    wsclean_fmem
        Maximum fraction of the total memory per node that WSClean may use (default = 0.9).

    ndir_per_node
        Maximum umber of directions to process in parallel on each node (default
        = 1). Note that the number of CPUs (set with the
        :term:`ncpu` parameter) and the amount of memory available to WSClean
        (set with the term:`wsclean_fmem` parameter) will be divided among the
        directions on each node.

.. _parset_checkfactor_options:

``[checkfactor]``
-----------------

.. glossary::

    facet_viewer
        Use ``casa`` or ``ds9`` for facet images (default = ``casa``).

    ds9_load_regions
        Load facet regions (ds9 only; default = ``False``).

    ds9_limits
        Scale limits (min max) in Jy/beam (ds9 only; default = full range).

    image_display
        Use ``display`` or ``eog`` to display PNG images (default = ``display``).


.. _parset_ms_specific_options:

``[<Your_MS_Name>]``
--------------------

MS-specific parameters (optional). You have to give the name of the MS (without
the path) as the section name. Currently, only the initial sky model can
be specified here.

.. glossary::

    init_skymodel
        Full path to the skymodel that was used to subtract the sources in the
	MS that was given as the section-name. For multi-epoch (interleaved or
	multi-night) observations the skymodel has to be specified only for one
	MS of each frequency group, it will then be used for all MSs in this
	frequency group. (Mixing MSs of the same frequency but in which different
	skymodels were used to subtract the sources is currently not possible.)
