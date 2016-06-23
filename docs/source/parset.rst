.. _factor_parset:

The Factor parset
=================

Before Factor can be run, a parset describing the reduction must be made. The
parset is a simple text file defining the parameters of a run in a number of
sections. For example, a typical parset for a basic reduction on a single
machine is shown below::

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
        1200). Generally, the number of chunks should be at least the number of
        available CPUs. To prevent Factor from chunking the data, set this value to be
        larger than the length of the longest dataset (in this case, Factor will not
        make copies of the files but will instead make symbolic links to them instead,
        so please make backup copies yourself).

    interactive
        Use interactive mode (default = ``False``). Factor will ask for confirmation of
        internally derived DDE calibrators and facets.

    keep_avg_facet_data
        Keep averaged calibrated data for each facet to allow re-imaging by hand
        (default = ``True``). If a target is specified (see below), the averaged
        data for the target is always kept, regardless of this setting. If the
        averaged data are kept, reimaging will be dramatically faster if
        multiple images per facet are made (e.g., at different scales)

    keep_unavg_facet_data
        Keep unaveraged calibrated data for each facet (default = ``False``).


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
        Size of frequency block in MHz over which a single TEC solution is fit
        (default = 10.0).

    peel_flux_Jy
        Peel the calibrator for sources above this flux density (default = 25.0).
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

    reimage_selfcaled
        Re-image all directions for which selfcal was successful (default = ``True``).

    wsclean_image_padding
        Padding factor for WSClean images (default = 1.6).

    wsclean_model_padding
        Padding factor for WSClean models (default = 1.4).

    max_peak_smearing
        Max desired peak flux density reduction at center of the facet edges due to
        bandwidth smearing (at the mean frequency) and time smearing (default = 0.15 =
        15% reduction in peak flux). Higher values result in shorter run times but
        more smearing away from the facet centers. This value only applies to the
        facet imaging (self calibration always uses a value of 0.15).

    facet_imager
        Use WSClean or CASA for imaging of entire facet (default = ``wsclean``). For large
        bandwidths, the CASA imager is typically faster.

    wsclean_nchannels_factor
        Max factor used to set the number of WSClean channel images when wide-band
        clean is used (default = 4). The number of channel images is determined by
        dividing the number of bands by the nearest divisor to this factor. Smaller
        values produce better results but require longer run times. Wide-band clean is
        activated when there are more than 5 bands.

    wsclean_add_bands
        Allow flagged data to be added during WSClean imaging to allow
        :term:`wsclean_nchannels_factor` to be a divisor of the number bands (default = ``True`).
        Enabling this option can dramatically speed up imaging with WSClean when the
        number of bands before padding does not allow :term:`wsclean_nchannels_factor` to be
        greater than 1 (e.g., :term:`wsclean_nchannels_factor` must be 1 to be an even divisor
        of 29 bands, so activating this option would add 1 band of flagged data to
        produce 30 bands, which will work with :term:`wsclean_nchannels_factor` = 3, 5, or 6)

    selfcal_cellsize_arcsec
        Self calibration pixel size in arcsec (default = 1.5).

    selfcal_robust
        Self calibration Briggs robust parameter for CASA (default = -0.25).

    selfcal_robust_wsclean
        Self calibration Briggs robust parameter for WSClean (default = -0.5).

    selfcal_min_uv_lambda
        Self calibration minimum uv distance in lambda (default = 80).

    selfcal_scales
        Self calibration multiscale clean scales (default = ``[0, 3, 7, 25, 60,
        150]``; set to ``[0]`` to disable multiscale clean).

    selfcal_clean_threshold
        Use a clean threshold during selfcal imaging (default = ``False``). If ``False``,
        clean will always stop at 1000 iterations. If ``True``, clean will stop when it
        reaches the 1 sigma noise level.

    selfcal_adaptive_threshold
        Use an adaptive masking threshold during selfcal imaging (default = ``False``). If
        ``True``, the masking threshold will be estimated using the negative peaks in the
        image, which can help selfcal convergence in the presence of strong artifacts.

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

    transfer_radius_deg
        Radius in degrees within which the direction-dependent solutions will be
        transferred before starting selfcal (default = 0.0; i.e., disabled). If a
        direction is within this distance of a calibrator for which selfcal was
        successful, the dir-dep selfcal solutions from this calibrator will be used
        instead of the dir-indep ones.

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
        PBS / torque reserved nodes, or use ``clusterdesc_file = JUROPA_slurm`` to use
        multiple nodes in a slurm reservation on JUROPA.
        If not given, the clusterdesc file for a single (i.e., local) node is used.

        .. note::

            On a cluster that uses torque and PBS, Factor will automatically determine the nodes for which you have a
            PBS reservation and use them. Note that you must ask for all the nodes you need
            in a single PBS script, so that all nodes are available for the full Factor run. An
            example PBS script is shown below::

                #!/bin/bash
                #PBS -N Factor
                #PBS -l walltime=100:00:00
                #PBS -l nodes=6:ppn=6

                cd $PBS_O_WORKDIR
                source ~rafferty/init_factor
                runfactor factor.parset

    dir_local
        Full path to a local disk on the nodes for I/O-intensive processing. The path
        must be the same for all nodes. Note: do not specify this parameter if you are
        running on a single machine, as it will cause conflicts between directions
        that are processed in parallel (no default).

    ncpu
        Maximum number of CPUs per node to use (default = all).

    wsclean_fmem
        Maximum fraction of the total memory per node that WSClean may use (default = 0.9).

    ndir_per_node
        Maximum umber of directions to process in parallel on each node (default
        = 1). Note that the number of CPUs (set with the
        :term:`ncpu` parameter) and the amount of memory available to WSClean
        :(set with the term:`wsclean_fmem` parameter) will be divided among the
        :directions on each node.


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
