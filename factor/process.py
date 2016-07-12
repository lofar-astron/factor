"""
Module that preforms the processing
"""
import sys
import os
import shutil
import numpy as np
import logging
import pickle
import collections
import casacore.tables as pt
from lofarpipe.support.data_map import DataMap
import factor
import factor.directions
import factor.parset
import factor.cluster
from factor.operations.outlier_ops import *
from factor.operations.field_ops import *
from factor.operations.facet_ops import *
from factor.lib.scheduler import Scheduler
from factor.lib.direction import Direction
from factor.lib.band import Band


log = logging.getLogger('factor')


def run(parset_file, logging_level='info', dry_run=False, test_run=False,
    reset_directions=[], reset_operations=[]):
    """
    Processes a dataset using facet calibration

    This function runs the operations in the correct order and handles all the
    bookkeeping for the processing: e.g., which operations are needed for each
    direction (depending on success of selfcal, the order in which they were
    processed, etc.).

    It also handles the set up of the computing parameters and the generation of
    DDE calibrators and facets.

    Parameters
    ----------
    parset_file : str
        Filename of parset containing processing parameters
    logging_level : str, optional
        One of 'degug', 'info', 'warning' in decreasing order of verbosity
    dry_run : bool, optional
        If True, do not run pipelines. All parsets, etc. are made as normal
    test_run : bool, optional
        If True, use test settings. These settings are for testing purposes
        only and will not produce useful results
    reset_directions : list of str, optional
        List of names of directions to be reset
    reset_operations : list of str, optional
        Llist of operations to be reset

    """
    # Read parset
    parset = factor.parset.parset_read(parset_file)

    # Set up logger
    parset['logging_level'] = logging_level
    factor._logging.set_level(logging_level)

    # Set up clusterdesc, node info, scheduler, etc.
    scheduler = _set_up_compute_parameters(parset, dry_run)

    # Prepare vis data
    bands = _set_up_bands(parset, test_run)

    # Set up directions and groups
    directions, direction_groups = _set_up_directions(parset, bands, dry_run,
    test_run, reset_directions, reset_operations)

    # Run peeling operations on outlier directions and any facets
    # for which the calibrator is to be peeled
    set_sub_data_colname = True
    peel_directions = [d for d in directions if d.is_outlier]
    peel_directions.extend([d for d in directions if d.peel_calibrator])
    if len(peel_directions) > 0:
        log.info('Peeling {0} direction(s)'.format(len(peel_directions)))

        # Do the peeling
        for d in peel_directions:
            # Reset if needed. Note that proper reset of the subtract steps in
            # outlierpeel and facetsub is not currently supported
            d.reset_state(['outlierpeel', 'facetpeel', 'facetpeelimage', 'facetsub'])

            if d.is_outlier:
                op = OutlierPeel(parset, bands, d)
            else:
                op = FacetPeel(parset, bands, d)
            scheduler.run(op)

            # Check whether direction went through selfcal successfully. If
            # not, exit
            if d.selfcal_ok:
                # Set the name of the subtracted data column for all directions
                if set_sub_data_colname:
                    for direction in directions:
                        direction.subtracted_data_colname = 'SUBTRACTED_DATA_ALL_NEW'
                    set_sub_data_colname = False
            else:
                log.error('Peeling failed for direction {0}.'.format(d.name))
                log.info('Exiting...')
                sys.exit(1)

            if d.peel_calibrator:
                # Do the imaging of the facet if calibrator was peeled and
                # subtract the improved model
                op = FacetPeelImage(parset, bands, d)
                scheduler.run(op)

                op = FacetSub(parset, bands, d)
                scheduler.run(op)

    # Run selfcal and subtract operations on direction groups
    for gindx, direction_group in enumerate(direction_groups):
        log.info('Self calibrating {0} direction(s) in Group {1}'.format(
            len(direction_group), gindx+1))

        # Set up reset of any directions that need it. If the direction has
        # already been through the facetsub operation, we must undo the
        # changes with the facetsubreset operation
        direction_group_reset = [d for d in direction_group if d.do_reset]
        direction_group_reset_facetsub = [d for d in direction_group_reset if
            'facetsub' in d.completed_operations]
        if (len(direction_group_reset_facetsub) > 0 and
            ('facetselfcal' in reset_operations or
            'facetsub' in reset_operations or
            'facetsubreset' in reset_operations)):
            for d in direction_group_reset_facetsub:
                if ('facetsubreset' in d.completed_operations or
                    'facetsubreset' in reset_operations):
                    # Reset a previous reset, but only if it completed successfully
                    # or is explicitly specified for reset (to allow one to resume
                    # facetsubreset instead of always resetting and restarting it)
                    d.reset_state('facetsubreset')
            ops = [FacetSubReset(parset, bands, d) for d in direction_group_reset_facetsub]
            for op in ops:
                scheduler.run(op)
        for d in direction_group_reset:
            d.reset_state(['facetselfcal', 'facetsub'])

        # Do selfcal or peeling on calibrator only
        ops = [FacetSelfcal(parset, bands, d) for d in direction_group]
        scheduler.run(ops)

        if dry_run:
            # For dryrun, skip selfcal verification
            for d in direction_group:
                d.selfcal_ok = True
        direction_group_ok = [d for d in direction_group if d.selfcal_ok]
        if set_sub_data_colname:
            # Set the name of the subtracted data column for remaining
            # directions (if needed)
            if len(direction_group_ok) > 0:
                for d in directions:
                    if d.name != direction_group_ok[0].name:
                        d.subtracted_data_colname = 'SUBTRACTED_DATA_ALL_NEW'
                set_sub_data_colname = False

        # Subtract final model(s) for directions for which selfcal went OK
        ops = [FacetSub(parset, bands, d) for d in direction_group_ok]
        for op in ops:
            scheduler.run(op)

        # Handle directions in this group for which selfcal failed
        selfcal_ok = [d.selfcal_ok for d in direction_group]
        for d in direction_group:
            if not d.selfcal_ok:
                log.warn('Self calibration failed for direction {0}.'.format(d.name))
        if not all(selfcal_ok) and parset['calibration_specific']['exit_on_selfcal_failure']:
            log.info('Exiting...')
            sys.exit(1)

    # Check that at least one direction went through selfcal successfully. If
    # not, exit
    if len([d for d in directions if d.selfcal_ok]) == 0:
        log.error('Self calibration failed for all directions. Exiting...')
        sys.exit(1)

    # (Re)image facets for each set of cellsize, robust, taper, and uv cut settings
    dirs_with_selfcal = [d for d in directions if d.selfcal_ok]
    dirs_with_selfcal_to_reimage = [d for d in dirs_with_selfcal if not d.is_patch
        and not d.is_outlier]
    dirs_without_selfcal = [d for d in directions if not d.selfcal_ok and not
        d.is_patch and not d.is_outlier]
    if len(dirs_without_selfcal) > 0:
        log.info('Imaging the following direction(s) with nearest self calibration solutions:')
        log.info('{0}'.format([d.name for d in dirs_without_selfcal]))
    for d in dirs_without_selfcal:
        # Search for nearest direction with successful selfcal
        nearest, sep = factor.directions.find_nearest(d, dirs_with_selfcal)
        log.debug('Using solutions from direction {0} for direction {1} '
            '(separation = {2} deg).'.format(nearest.name, d.name, sep))
        d.dir_dep_parmdb_mapfile = nearest.dir_dep_parmdb_mapfile
        d.save_state()

    if len(dirs_with_selfcal_to_reimage + dirs_without_selfcal) > 0:
        cellsizes = parset['imaging_specific']['facet_cellsize_arcsec']
        tapers = parset['imaging_specific']['facet_taper_arcsec']
        robusts = parset['imaging_specific']['facet_robust']
        min_uvs = parset['imaging_specific']['facet_min_uv_lambda']
        nimages = len(cellsizes)

        for image_indx, (cellsize_arcsec, taper_arcsec, robust, min_uv_lambda) in enumerate(
            zip(cellsizes, tapers, robusts, min_uvs)):
            # Always image directions that did not go through selfcal
            dirs_to_image = dirs_without_selfcal[:]

            # Only reimage facets with selfcal imaging parameters if reimage_selfcal flag is set
            if parset['imaging_specific']['facet_imager'] == 'wsclean':
                selfcal_robust = parset['imaging_specific']['selfcal_robust_wsclean']
            else:
                selfcal_robust = parset['imaging_specific']['selfcal_robust']
            if (cellsize_arcsec == parset['imaging_specific']['selfcal_cellsize_arcsec'] and
                robust == selfcal_robust and taper_arcsec == 0.0):
                if parset['imaging_specific']['reimage_selfcaled']:
                    dirs_to_image += dirs_with_selfcal_to_reimage
            else:
                dirs_to_image += dirs_with_selfcal_to_reimage

            if len(dirs_to_image) > 0:
                log.info('Imaging with cellsize = {0} arcsec, robust = {1}, '
                    'taper = {2} arcsec, min_uv = {3} lambda'.format(cellsize_arcsec,
                    robust, taper_arcsec, min_uv_lambda))
                log.info('Imaging the following direction(s):')
                log.info('{0}'.format([d.name for d in dirs_to_image]))

            # Reset facetimage op for any directions that need it
            directions_reset = [d for d in dirs_to_image if d.do_reset]
            for d in directions_reset:
                op = FacetImage(parset, bands, d, cellsize_arcsec, robust,
                    taper_arcsec, min_uv_lambda)
                d.reset_state(op.name)

            # Do facet imaging
            ops = [FacetImage(parset, bands, d, cellsize_arcsec, robust,
                taper_arcsec, min_uv_lambda) for d in dirs_to_image]
            scheduler.run(ops)

        # Mosaic the final facet images together
        for i, (cellsize_arcsec, taper_arcsec, robust, min_uv_lambda) in enumerate(
            zip(cellsizes, tapers, robusts, min_uvs)):
            if parset['imaging_specific']['make_mosaic']:
                # Make direction object for the field and load previous state (if any)
                field = Direction('field', bands[0].ra, bands[0].dec,
                    factor_working_dir=parset['dir_working'])
                field.load_state()

                # Set averaging for primary beam generation
                field.avgpb_freqstep = bands[0].nchan
                field.avgpb_timestep = int(120.0 / bands[0].timepersample)

                # Reset the field direction if specified
                if 'field' in reset_directions:
                    op = FieldMosaic(parset, bands, field, cellsize_arcsec, robust,
                        taper_arcsec, min_uv_lambda)
                    field.reset_state(op.name)

                field.facet_image_filenames = []
                field.facet_vertices_filenames = []
                for d in directions:
                    if not d.is_patch:
                        facet_image = DataMap.load(d.facet_image_mapfile)[0].file
                        field.facet_image_filenames.append(facet_image)
                        field.facet_vertices_filenames.append(d.save_file)

                # Do mosaicking
                op = FieldMosaic(parset, bands, field, cellsize_arcsec, robust,
                        taper_arcsec, min_uv_lambda)
                scheduler.run(op)

    log.info("Factor has finished :)")


def _set_up_compute_parameters(parset, dry_run=False):
    """
    Sets up compute parameters and operation scheduler

    Parameters
    ----------
    parset : dict
        Parset containing processing parameters
    dry_run : bool, optional
        If True, do not run pipelines. All parsets, etc. are made as normal

    Returns
    -------
    scheduler : Scheduler instance
        The operation scheduler used by the run() function

    """
    log.info('Setting up cluster/node parameters...')

    cluster_parset = parset['cluster_specific']
    if 'clusterdesc_file' not in cluster_parset:
        log.warn('Did not find \"clusterdesc_file\" in parset-dict! This shouldn\'t happen, but trying to continue anyhow.')
        parset['cluster_specific']['clusterdesc'] = 'local.clusterdesc'
    else:
        if (cluster_parset['clusterdesc_file'].lower() == 'pbs' or
            ('cluster_type' in cluster_parset and cluster_parset['cluster_type'].lower() == 'pbs')):
            log.info('Using cluster setting: "PBS".')
            parset['cluster_specific']['clusterdesc'] = factor.cluster.make_pbs_clusterdesc()
            parset['cluster_specific']['clustertype'] = 'pbs'
        elif (cluster_parset['clusterdesc_file'].lower() == 'juropa_slurm' or
            ('cluster_type' in cluster_parset and cluster_parset['cluster_type'].lower() == 'juropa_slurm')):
            log.info('Using cluster setting: "JUROPA_slurm" (Single genericpipeline using multiple nodes).')
            # slurm_srun on JUROPA uses the local.clusterdesc
            parset['cluster_specific']['clusterdesc'] = os.path.join(parset['lofarroot'], 'share', 'local.clusterdesc')
            parset['cluster_specific']['clustertype'] = 'juropa_slurm'
            parset['cluster_specific']['node_list'] = ['localhost']
        else:
            log.info('Using cluster setting: "local" (Single node).')
            parset['cluster_specific']['clusterdesc'] = cluster_parset['clusterdesc_file']
            parset['cluster_specific']['clustertype'] = 'local'
    if not 'node_list' in parset['cluster_specific']:
        parset['cluster_specific']['node_list'] = factor.cluster.get_compute_nodes(
            parset['cluster_specific']['clusterdesc'])

    # check ulimit(s)
    try:
        import resource
        nof_files_limits = resource.getrlimit(resource.RLIMIT_NOFILE)
        if parset['cluster_specific']['clustertype'] == 'local' and nof_files_limits[0] < nof_files_limits[1]:
            log.debug('Setting limit for number of open files to: {}.'.format(nof_files_limits[1]))
            resource.setrlimit(resource.RLIMIT_NOFILE,(nof_files_limits[1],nof_files_limits[1]))
            nof_files_limits = resource.getrlimit(resource.RLIMIT_NOFILE)
        log.debug('Active limit for number of open files is {0}, maximum limit is {1}.'.format(nof_files_limits[0],nof_files_limits[1]))
        if nof_files_limits[0] < 2048:
            log.warn('The limit for number of open files is small, this could results in a "Too many open files" problem when running factor.')
            log.warn('The active limit can be increased to the maximum for the user with: "ulimit -Sn <number>" (bash) or "limit descriptors 1024" (csh).')
    except resource.error:
        log.warn('Cannot check limits for number of open files, what kind of system is this?')

    # Get paths to required executables
    factor.cluster.find_executables(parset)

    # Set up scheduler for operations (pipeline runs)
    ndir_simul = len(parset['cluster_specific']['node_list']) * \
        parset['cluster_specific']['ndir_per_node']
    if parset['direction_specific']['groupings'] is not None:
        ngroup_max = int(max([int(n.items()[0][0]) for n in
            parset['direction_specific']['groupings']]))
    else:
        ngroup_max = 1
    if ndir_simul < ngroup_max:
        log.warn('The maximum number of directions that can be proccessed '
            'simultaneously ({0}) is less than the number of directions in the '
            'largest group ({1}). For best performance, these values should be '
            'equal'.format(ndir_simul, ngroup_max))
    scheduler = Scheduler(parset['genericpipeline_executable'], max_procs=ndir_simul,
        dry_run=dry_run)

    return scheduler


def _set_up_bands(parset, test_run=False):
    """
    Sets up bands for processing

    Parameters
    ----------
    parset : dict
        Parset containing processing parameters
    test_run : bool, optional
        If True, use test settings. These settings are for testing purposes
        only and will not produce useful results

    Returns
    -------
    bands : List of Band instances
        All bands to be used by the run() function, ordered from low to high
        frequency

    """
    log.info('Checking input bands...')
    msdict = {}
    for ms in parset['mss']:
        # group all found MSs by frequency
        sw = pt.table(ms+'::SPECTRAL_WINDOW', ack=False)
        msfreq = int(sw.col('REF_FREQUENCY')[0])
        sw.close()
        if msfreq in msdict:
            msdict[msfreq].append(ms)
        else:
            msdict[msfreq] = [ms]
    bands = []
    for MSkey in msdict.keys():
        # Check for any sky models specified by user
        # there only needs to be a skymodel specyfied for one file in each band
        skymodel_dirindep = None
        for ms in msdict[MSkey]:
            msbase = os.path.basename(ms)
            if msbase in parset['ms_specific']:
                if 'init_skymodel' in parset['ms_specific'][msbase]:
                    skymodel_dirindep = parset['ms_specific'][msbase]['init_skymodel']
                    if not os.path.exists(skymodel_dirindep):
                        log.error('Sky model specified in parset for band {} was '
                            'not found. Exiting...'.format(msbase))
                        sys.exit(1)
                    break
        band = Band(msdict[MSkey], parset['dir_working'], parset['parmdb_name'],
            skymodel_dirindep, local_dir=parset['cluster_specific']['dir_local'],
            test_run=test_run, chunk_size_sec=parset['chunk_size_sec'])
        bands.append(band)

    # Sort bands by frequency
    band_freqs = [band.freq for band in bands]
    bands = np.array(bands)[np.argsort(band_freqs)].tolist()

    # Check bands for problems
    nchan_list = []
    ra_list = []
    dec_list = []
    has_gaps = False
    for band in bands:
        if len(band.missing_channels) > 0:
            log.error('Found one or more frequency gaps in band {}'.format(band.msnames[0]))
            has_gaps = True
        nchan_list.append(band.nchan)
        ra_list.append(band.ra)
        dec_list.append(band.dec)

    # Check for frequency gaps
    if has_gaps:
        log.error('Bands cannot have frequency gaps. Exiting...')
        sys.exit(1)

    # Check that all bands have the same number of channels
    duplicate_chans = set(nchan_list)
    if len(duplicate_chans) != 1:
        for d in duplicate_chans:
            bands_with_duplicates = [band.msnames[0] for band in bands if band.nchan == d]
            log.error('Found {0} channels in band(s): {1}'.format(d,
                ', '.join(bands_with_duplicates)))
        log.error('All bands must have the same number of channels. Exiting...')
        sys.exit(1)

    # Check that number of channels supports enough averaging steps
    nchans = list(duplicate_chans)[0]
    n_divisors = len([x+1 for x in range(nchans) if not nchans % (x+1)])
    if n_divisors < 5:
        log.warn('Number of channels per band is {}. Averaging will '
            'not work well (too few divisors)'.format(nchans))

    # Check that phase center is the same for all bands
    pos_tol_deg = 1.0 / 3600.0 # positional tolerance
    for ra, dec in zip(ra_list[1:], dec_list[1:]):
        if (not factor.directions.approx_equal(ra, ra_list[0], tol=pos_tol_deg) or
            not factor.directions.approx_equal(dec, dec_list[0], tol=pos_tol_deg)):
            log.error('Input bands do not have a common phase center. Exiting...')
            sys.exit(1)

    # Determine whether any bands lack an initial sky model
    bands_no_skymodel = [b for b in bands if b.skymodel_dirindep is None]
    if len(bands_no_skymodel) > 0:
        if len(bands_no_skymodel) > 0:
            band_names = ', '.join([b.name for b in bands_no_skymodel])
            log.error('A direction-indpendent sky model was not found for the '
                'following bands: {}'.format(band_names))
        log.info('Exiting...')
        sys.exit(1)

    return bands


def _set_up_directions(parset, bands, dry_run=False, test_run=False,
    reset_directions=[], reset_operations=[]):
    """
    Sets up directions (facets)

    Parameters
    ----------
    parset : dict
        Parset containing processing parameters
    bands : list of Band instances
        Vis data
    dry_run : bool, optional
        If True, do not run pipelines. All parsets, etc. are made as normal
    test_run : bool, optional
        If True, use test settings. These settings are for testing purposes
        only and will not produce useful results
    reset_directions : list of str, optional
        List of direction names to be reset
    reset_operations : list of str, optional
        Llist of operations to be reset

    Returns
    -------
    directions : List of Direction instances
        All directions to be used by the run() function
    direction_groups : List of lists of Direction instances
        Groups of directions to be selfcal-ed

    """
    dir_parset = parset['direction_specific']

    max_radius_deg = dir_parset['max_radius_deg']
    ref_band = bands[-1]

    if dir_parset['faceting_skymodel'] is not None:
        import lsmtool
        log.info("Using {} as sky model for source avoidance and DDE calibrator "
            "selection (if desired)".format(dir_parset['faceting_skymodel']))
        initial_skymodel = lsmtool.load(dir_parset['faceting_skymodel'])
    else:
        log.info("Building local sky model for source avoidance and DDE calibrator "
            "selection (if desired)...")
        initial_skymodel = factor.directions.make_initial_skymodel(ref_band)

    log.info('Setting up directions...')
    directions = _initialize_directions(parset, initial_skymodel, ref_band,
        max_radius_deg=max_radius_deg, dry_run=dry_run)

    # Check with user
    if parset['interactive']:
        print("\nFacet and DDE calibrator regions saved. Please check that they\n"
            "are OK before continuing. You can edit the directions file and\n"
            "continue; FACTOR will pick up any changes to it. Note: if you\n"
            "choose not to continue and you let FACTOR generate the directions\n"
            "internally, you must delete the FACTOR-made directions file\n"
            "(dir_working/factor_directions.txt) before restarting if you want\n"
            "to FACTOR to regenerate it\n")
        prompt = "Continue processing (y/n)? "
        answ = raw_input(prompt)
        while answ.lower() not in  ['y', 'n', 'yes', 'no']:
            answ = raw_input(prompt)
        if answ.lower() in ['n', 'no']:
            log.info('Exiting...')
            sys.exit(0)
        else:
            # Continue processing, but first re-initialize the directions to
            # pick up any changes the user made to the directions file
            directions = _initialize_directions(parset, initial_skymodel,
                ref_band, max_radius_deg=max_radius_deg, dry_run=dry_run)

    # Warn user if they've specified a direction to reset that does not exist
    direction_names = [d.name for d in directions]
    for name in reset_directions:
        if name not in direction_names and name != 'field':
            log.warn('Direction {} was specified for resetting but does not '
                'exist in current list of directions'.format(name))

    # Load previously completed operations (if any) and facetsubreset-specific
    # attributes and save the state
    for direction in directions:
        direction.load_state()
        direction.save_state()

    # Select subset of directions to process
    target_has_own_facet = dir_parset['target_has_own_facet']
    if target_has_own_facet:
        direction_names = [d.name for d in directions]
        target = directions[direction_names.index('target')]
    if dir_parset['ndir_process'] is not None:
        if dir_parset['ndir_process'] < len(directions):
            directions = directions[:dir_parset['ndir_process']]

            # Make sure target is still included
            direction_names = [d.name for d in directions]
            if target_has_own_facet and 'target' not in direction_names:
                directions.append(target)

            # Warn user if reimaging is to be done but they are not processing
            # the full field
            if parset['imaging_specific']['reimage_selfcaled']:
                log.warn("The reimage_selfcaled parameter is True but all directions "
                    "will not be processed. If you're interested in only a single "
                    "target in the last facet, then re-imaging will not improve results.")

    # Set various direction attributes
    for i, direction in enumerate(directions):
        # Set direction sky model
        direction.set_skymodel(initial_skymodel.copy())

        # Set peeling flag (i.e., facet calibrator should be peeled before facet
        # is imaged)
        total_flux_jy, peak_flux_jy_bm = direction.get_cal_fluxes()
        effective_flux_jy = peak_flux_jy_bm * (total_flux_jy / peak_flux_jy_bm)**0.667
        if (effective_flux_jy > parset['calibration_specific']['peel_flux_jy'] or
            direction.is_outlier):
            direction.find_peel_skymodel()
            if direction.peel_skymodel is not None:
                if not direction.is_outlier:
                   direction.peel_calibrator = True
                log.info('Direction {0} will be peeled using sky model: {1}'.format(
                    direction.name, direction.peel_skymodel))
            else:
                if direction.is_outlier:
                    log.error('Direction {} was specified as an outlier source '
                    'but an appropriate sky model is not available'.format(direction.name))
                    sys.exit(1)
                else:
                    log.warning('The flux density of direction {} exceeds peel_flux_Jy '
                        'but an appropriate peeling sky model is not available. '
                        'This direction will go through normal self calibration '
                        'instead'.format(direction.name))

        # Set full correlation solve
        if effective_flux_jy > parset['calibration_specific']['solve_all_correlations_flux_jy']:
            if not parset['calibration_specific']['spline_smooth2d']:
                log.error('The option spline_smooth2d must be enabled to use '
                    'XY and YX correlations during the slow gain solve')
                sys.exit(1)
            direction.solve_all_correlations = True

        # Set skip_facet_imaging flag
        direction.skip_facet_imaging = parset['imaging_specific']['skip_facet_imaging']

        # Set field center to that of first band (all bands have the same phase
        # center)
        direction.field_ra = bands[0].ra
        direction.field_dec = bands[0].dec

        # Reset state if specified
        if direction.name in reset_directions:
            direction.do_reset = True
            if len(reset_operations) > 0:
                direction.reset_operations = reset_operations
            else:
                direction.reset_operations = (direction.completed_operations[:] +
                    direction.started_operations[:])
        else:
            direction.do_reset = False

    # Select directions to selfcal, excluding outliers and target
    if target_has_own_facet:
        # Make sure target is not a DDE calibrator and is at end of directions list
        selfcal_directions = [d for d in directions if d.name != target.name and
                              not d.is_outlier and not d.peel_calibrator]
        directions = [d for d in directions if d.name != target.name] + [target]
    else:
        selfcal_directions = [d for d in directions if not d.is_outlier and
            not d.peel_calibrator]

    if dir_parset['ndir_selfcal'] is not None:
        if dir_parset['ndir_selfcal'] <= len(selfcal_directions):
            selfcal_directions = selfcal_directions[:dir_parset['ndir_selfcal']]

    # Divide directions into groups for selfcal
    direction_groups = factor.directions.group_directions(selfcal_directions,
        n_per_grouping=dir_parset['groupings'], allow_reordering=dir_parset['allow_reordering'])

    return directions, direction_groups


def _initialize_directions(parset, initial_skymodel, ref_band, max_radius_deg=None,
    dry_run=False):
    """
    Read in directions file and initialize resulting directions

    Parameters
    ----------
    parset : dict
        Parset containing processing parameters
    initial_skymodel : SkyModel object
        Local sky model
    ref_band : Band object
        Reference band
    max_radius_deg : float, optional
        Maximum radius in degrees from the phase center within which to include
        sources. If None, it is set to the FWHM (i.e., a diameter of 2 * FWHM)
    dry_run : bool, optional
        If True, do not run pipelines. All parsets, etc. are made as normal

    Returns
    -------
    directions : List of Direction instances
        All directions to be used

    """
    dir_parset = parset['direction_specific']
    s = initial_skymodel.copy()

    # First check for user-supplied directions file, then for Factor-generated
    # file from a previous run, then for parameters needed to generate it internally
    if 'directions_file' in dir_parset:
        directions = factor.directions.directions_read(dir_parset['directions_file'],
            parset['dir_working'])
    elif os.path.exists(os.path.join(parset['dir_working'], 'factor_directions.txt')):
        directions = factor.directions.directions_read(os.path.join(parset['dir_working'],
            'factor_directions.txt'), parset['dir_working'])
    else:
        if dir_parset['flux_min_jy'] is None or \
            dir_parset['size_max_arcmin'] is None or \
            dir_parset['separation_max_arcmin'] is None:
                log.critical('If no directions file is specified, you must '
                    'give values for flux_min_Jy, size_max_arcmin, and '
                    'separation_max_arcmin')
                sys.exit(1)
        else:
            # Make directions from dir-indep sky model of highest-frequency
            # band, as it has the smallest field of view
            log.info("No directions file given. Selecting directions internally...")

            # Filter out sources that lie outside of maximum specific radius from phase
            # center
            if max_radius_deg is None:
                max_radius_deg = ref_band.fwhm_deg # means a diameter of 2 * FWHM
            log.info('Removing sources beyond a radius of {0} degrees (corresponding to '
                'a diameter of {1} * FWHM of the primary beam at {2} MHz)...'.format(
                max_radius_deg, round(2.0*max_radius_deg/ref_band.fwhm_deg, 1), ref_band.freq/1e6))
            dist = s.getDistance(ref_band.ra, ref_band.dec, byPatch=True)
            s.remove(dist > max_radius_deg, aggregate=True)

            # Generate the directions file
            if dir_parset['minimize_nonuniformity']:
                dir_parset['directions_file'] = factor.directions.make_directions_file_from_skymodel_uniform(
                    s, dir_parset['flux_min_jy'], dir_parset['size_max_arcmin'],
                    dir_parset['separation_max_arcmin'],
                    directions_max_num=dir_parset['ndir_max'],
                    interactive=parset['interactive'], ncpu=parset['cluster_specific']['ncpu'],
                    flux_min_for_merging_Jy=dir_parset['flux_min_for_merging_jy'])
            else:
                dir_parset['directions_file'] = factor.directions.make_directions_file_from_skymodel(
                    s, dir_parset['flux_min_jy'], dir_parset['size_max_arcmin'],
                    dir_parset['separation_max_arcmin'],
                    directions_max_num=dir_parset['ndir_max'],
                    interactive=parset['interactive'],
                    flux_min_for_merging_Jy=dir_parset['flux_min_for_merging_jy'])
            directions = factor.directions.directions_read(dir_parset['directions_file'],
                parset['dir_working'])

    # Add the target to the directions list if desired
    target_ra = dir_parset['target_ra']
    target_dec = dir_parset['target_dec']
    target_radius_arcmin = dir_parset['target_radius_arcmin']
    target_has_own_facet = dir_parset['target_has_own_facet']
    if target_has_own_facet:
        if target_ra is not None and target_dec is not None and target_radius_arcmin is not None:
            # Make target object
            target = Direction('target', target_ra, target_dec,
                factor_working_dir=parset['dir_working'])

            # Check if target is already in directions list because it was
            # selected as a DDE calibrator. If so, remove the duplicate
            nearest, dist = factor.directions.find_nearest(target, directions)
            if dist < dir_parset['target_radius_arcmin']/60.0:
                directions.remove(nearest)

            # Add target to directions list
            directions.append(target)
        else:
            log.critical('target_has_own_facet = True, but target RA, Dec, or radius not found in parset')
            sys.exit(1)

    # Set calibrator size (must be done before faceting below is done)
    for d in directions:
        d.set_cal_size(parset['imaging_specific']['selfcal_cellsize_arcsec'])

    # Create facets and patches
    faceting_radius_deg = dir_parset['faceting_radius_deg']
    if faceting_radius_deg is None:
        faceting_radius_deg = 1.25 * ref_band.fwhm_deg / 2.0
    beam_ratio = 1.0 / np.sin(ref_band.mean_el_rad) # ratio of N-S to E-W beam
    factor.directions.thiessen(directions, ref_band.ra, ref_band.dec,
        faceting_radius_deg, s=s, check_edges=dir_parset['check_edges'],
        target_ra=target_ra, target_dec=target_dec,
        target_radius_arcmin=target_radius_arcmin, beam_ratio=beam_ratio)

    # Make DS9 region files so user can check the facets, etc.
    ds9_facet_reg_file = os.path.join(parset['dir_working'], 'regions', 'facets_ds9.reg')
    factor.directions.make_ds9_region_file(directions, ds9_facet_reg_file)
    ds9_calimage_reg_file = os.path.join(parset['dir_working'], 'regions', 'calimages_ds9.reg')
    factor.directions.make_ds9_calimage_file(directions, ds9_calimage_reg_file)

    return directions
