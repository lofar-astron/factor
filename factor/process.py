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
from lofarpipe.support.data_map import DataMap
import factor
import factor.directions
import factor.parset
import factor.cluster
from factor.operations.field_ops import *
from factor.operations.facet_ops import *
from factor.lib.scheduler import Scheduler
from factor.lib.direction import Direction
from factor.lib.band import Band


log = logging.getLogger('factor')


def run(parset_file, logging_level='info', dry_run=False, test_run=False,
    reset_directions=[]):
    """
    Processes a dataset using facet calibration

    This function runs the operations in the correct order and handles all the
    bookkeeping for the processing: e.g., which operations are needed for each
    direction (depending on sucess of selfcal, the order in which they were
    processed, etc.).

    It also handles the setup of the computing parameters and the generation of
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
        List of direction names to be reset

    """
    # Read parset
    parset = factor.parset.parset_read(parset_file)

    # Set up logger
    parset['logging_level'] = logging_level
    factor._logging.set_level(logging_level)

    # Set up clusterdesc, node info, scheduler, etc.
    scheduler = _set_up_compute_parameters(parset, dry_run)

    # Prepare vis data
    bands, bands_initsubtract = _set_up_bands(parset, test_run)

    # Make direction object for the field
    field = Direction('field', bands[0].ra, bands[0].dec,
        factor_working_dir=parset['dir_working'])
    field.set_averaging_steps_and_solution_intervals(bands[0].chan_width_hz, bands[0].nchan,
        bands[0].timepersample, bands[0].nsamples, len(bands), parset['preaverage_flux_jy'])

    # Run initial sky model generation and create empty datasets
    if len(bands_initsubtract) > 0:
        # Reset the field direction if specified
        if 'field' in reset_directions:
            field.reset_state('initsubtract')

        log.info('Running initsubtract operation for bands: {0}'.
            format([b.name for b in bands_initsubtract]))
        # Combine the nodes and cores for the serial subtract operations
        field = factor.cluster.combine_nodes([field],
            parset['cluster_specific']['node_list'],
            parset['cluster_specific']['nimg_per_node'],
            parset['cluster_specific']['ncpu'],
            parset['cluster_specific']['fmem'],
            len(bands_initsubtract))[0]
        op = InitSubtract(parset, bands_initsubtract, field)
        scheduler.run(op)
    else:
        log.info("Sky models and SUBTRACTED_DATA_ALL found for all bands. "
            "Skipping initsubtract operation...")

    # Remove bands that failed initsubtract operation
    bands = [b for b in bands if not b.skip]

    # Define directions
    directions, direction_groups = _set_up_directions(parset, bands, field,
        dry_run, test_run, reset_directions)

    # Run selfcal and subtract operations on directions
    set_sub_data_colname = True
    for gindx, direction_group in enumerate(direction_groups):
        log.info('Processing {0} direction(s) in Group {1}'.format(
            len(direction_group), gindx+1))

        # Set up reset of any directions that need it. If the direction has
        # already been through the facetsub operation, we must undo the
        # changes with the facetsubreset operation before we reset facetselfcal
        # (otherwise the model data required to reset facetsub will be deleted)
        direction_group_reset = [d for d in direction_group if d.do_reset]
        direction_group_reset_facetsub = [d for d in direction_group_reset if
            'facetsub' in d.completed_operations]
        if len(direction_group_reset_facetsub) > 0:
            for d in direction_group_reset_facetsub:
                d.reset_state('facetsubreset')
            direction_group_reset_facetsub = factor.cluster.combine_nodes(
                direction_group_reset_facetsub,
                parset['cluster_specific']['node_list'],
                parset['cluster_specific']['nimg_per_node'],
                parset['cluster_specific']['ncpu'],
                parset['cluster_specific']['fmem'],
                len(bands))
            ops = [FacetSubReset(parset, bands, d) for d in direction_group_reset_facetsub]
            for op in ops:
                scheduler.run(op)
        for d in direction_group_reset:
            d.reset_state(['facetselfcal', 'facetsub'])

        # Divide up the nodes and cores among the directions for the parallel
        # selfcal operations
        direction_group = factor.cluster.divide_nodes(direction_group,
            parset['cluster_specific']['node_list'],
            parset['cluster_specific']['ndir_per_node'],
            parset['cluster_specific']['nimg_per_node'],
            parset['cluster_specific']['ncpu'],
            parset['cluster_specific']['fmem'],
            len(bands))

        # Check for any directions within transfer radius that have successfully
        # gone through selfcal
        dirs_with_selfcal = [d for d in directions if d.selfcal_ok]
        if len(dirs_with_selfcal) > 0:
            for d in direction_group:
                nearest, sep = factor.directions.find_nearest(d, dirs_with_selfcal)
                if sep < parset['direction_specific']['transfer_radius']:
                    log.debug('Initializing selfcal for direction {0} with solutions from direction {1}.'.format(
                        d.name, nearest.name))
                    d.dir_dep_parmdb_mapfile = nearest.dir_dep_parmdb_mapfile
                    d.save_state()
                    d.transfer_nearest_solutions = True

        # Do selfcal on calibrator only
        ops = [FacetSelfcal(parset, bands, d) for d in direction_group]
        scheduler.run(ops)
        if dry_run:
            # For dryrun, skip check
            for d in direction_group:
                d.selfcal_ok = True
        direction_group_ok = [d for d in direction_group if d.selfcal_ok]
        if set_sub_data_colname:
            if len(direction_group_ok) > 0:
                for d in directions:
                    if d.name != direction_group_ok[0].name:
                        d.subtracted_data_colname = 'SUBTRACTED_DATA_ALL_NEW'
                set_sub_data_colname = False

        # Combine the nodes and cores for the serial subtract operations
        direction_group_ok = factor.cluster.combine_nodes(direction_group_ok,
            parset['cluster_specific']['node_list'],
            parset['cluster_specific']['nimg_per_node'],
            parset['cluster_specific']['ncpu'],
            parset['cluster_specific']['fmem'],
            len(bands))

        # Subtract final model(s) for directions for which selfcal went OK
        ops = [FacetSub(parset, bands, d) for d in direction_group_ok]
        for op in ops:
            scheduler.run(op)

        # Handle directions in this group for which selfcal failed
        selfcal_ok = [d.selfcal_ok for d in direction_group]
        for d in direction_group:
            if not d.selfcal_ok:
                log.warn('Selfcal verification failed for direction {0}.'.format(d.name))
        if not all(selfcal_ok) and parset['exit_on_selfcal_failure']:
            log.info('Exiting...')
            sys.exit(1)

    # Check that at least one direction went through selfcal successfully. If
    # not, exit
    if len([d for d in directions if d.selfcal_ok]) == 0:
        log.error('Selfcal verification failed for all directions. Exiting...')
        sys.exit(1)

    # Make final facet images (from final empty datasets) if desired. Also image
    # any facets for which selfcal failed or no selfcal was done
    #
    # TODO: combine new facet sky models and adjust facet edges for new sources
    # (but only if all facets are to be re-imaged)
    #
    dirs_to_image = [d for d in directions if d.make_final_image and d.selfcal_ok]
    if len(dirs_to_image) > 0:
        log.info('Reimaging the following direction(s):')
        log.info('{0}'.format([d.name for d in dirs_to_image]))

    # Add directions without selfcal
    dirs_to_transfer = [d for d in directions if not d.selfcal_ok]
    if len(dirs_to_transfer) > 0:
        log.info('Imaging the following direction(s) with nearest selcal solutions:')
        log.info('{0}'.format([d.name for d in dirs_to_transfer]))
    dirs_with_selfcal = [d for d in directions if d.selfcal_ok]
    for d in dirs_to_transfer:
        # Search for nearest direction with successful selfcal
        nearest, sep = factor.directions.find_nearest(d, dirs_with_selfcal)
        log.debug('Using solutions from direction {0} for direction {1}.'.format(
            nearest.name, d.name))
        d.dir_dep_parmdb_mapfile = nearest.dir_dep_parmdb_mapfile
        d.save_state()
    dirs_to_image.extend(dirs_to_transfer)

    if len(dirs_to_image) > 0:
        # Set up reset of any directions that need it
        directions_reset = [d for d in dirs_to_image if d.do_reset]
        for d in directions_reset:
            d.reset_state('facetimage')

        # Group directions. This is done to ensure that multiple directions
        # aren't competing for the same resources
        ndir_simul = (len(parset['cluster_specific']['node_list']) *
            parset['cluster_specific']['ndir_per_node'])
        for i in range(int(np.ceil(len(dirs_to_image)/float(ndir_simul)))):
            dir_group = dirs_to_image[i*ndir_simul:(i+1)*ndir_simul]

            # Divide up the nodes and cores among the directions for the parallel
            # imaging operations
            dir_group = factor.cluster.divide_nodes(dir_group,
                parset['cluster_specific']['node_list'],
                parset['cluster_specific']['ndir_per_node'],
                parset['cluster_specific']['nimg_per_node'],
                parset['cluster_specific']['ncpu'],
                parset['cluster_specific']['fmem'],
                len(bands))

            ops = [FacetImage(parset, bands, d) for d in dir_group]
            scheduler.run(ops)

    # Mosaic the final facet images together
    if parset['make_mosaic']:
        # Reset the field direction if specified
        if 'field' in reset_directions or 'mosaic' in reset_directions:
            field.reset_state('makemosaic')

        field.facet_image_filenames = []
        field.facet_vertices_filenames = []
        for d in directions:
            facet_image = DataMap.load(d.facet_image_mapfile)[0].file
            field.facet_image_filenames.append(facet_image)
            field.facet_vertices_filenames.append(d.save_file)
        op = MakeMosaic(parset, bands, field)
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
        parset['cluster_specific']['clusterdesc'] = 'local.clusterdesc'
    else:
        if cluster_parset['clusterdesc_file'].lower() == 'pbs':
            parset['cluster_specific']['clusterdesc'] = factor.cluster.make_pbs_clusterdesc()
        else:
            parset['cluster_specific']['clusterdesc'] = cluster_parset['clusterdesc_file']
    if not 'node_list' in cluster_parset:
        parset['cluster_specific']['node_list'] = factor.cluster.get_compute_nodes(
            parset['cluster_specific']['clusterdesc'])

    # Get paths to required executables
    factor.cluster.find_executables(parset)

    # Set up scheduler for operations (pipeline runs)
    ndir_simul = len(parset['cluster_specific']['node_list']) * \
        parset['cluster_specific']['ndir_per_node']
    if parset['direction_specific']['groupings'] is not None:
        ngroup_max = int(max(parset['direction_specific']['groupings'].keys()))
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
    bands_initsubtract : List of Band instances
        Subset of bands for InitSubtract operation, ordered from low to high
        frequency

    """
    log.info('Checking input bands...')
    bands = []
    for ms in parset['mss']:
        # Check for any sky models specified by user
        skymodel_dirindep = None
        msbase = os.path.basename(ms)
        if msbase in parset['ms_specific']:
            if 'init_skymodel' in parset['ms_specific'][msbase]:
                skymodel_dirindep = parset['ms_specific'][msbase]['init_skymodel']
                if not os.path.exists(skymodel_dirindep):
                    log.error('Sky model specified in parset for band {} was '
                        'not found. Exiting...'.format(msbase))
                    sys.exit(1)

        band = Band(ms, parset['dir_working'], parset['parmdb_name'], skymodel_dirindep,
            test_run=test_run)
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
            log.error('Found one or more frequency gaps in band {}'.format(band.msname))
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
            bands_with_duplicates = [band.msname for band in bands if band.nchan == d]
            log.error('Found {0} channels in band(s): {1}'.format(d,
                ', '.join(bands_with_duplicates)))
        log.error('All bands must have the same number of channels. Exiting...')
        sys.exit(1)

    # Check that number of channels supports enough averaging steps
    if list(duplicate_chans)[0] not in [18, 20, 24]:
        log.warn('Number of channels per band is not 18, 20, or 24. Averaging will '
            'not work well (too few divisors)')

    # Determine whether any bands need to be run through the initsubract operation.
    # This operation is only needed if band lacks an initial skymodel or
    # the SUBTRACTED_DATA_ALL column
    bands_initsubtract = []
    for band in bands:
        if band.skymodel_dirindep is None or not band.has_sub_data:
            # Set image sizes for initsubtract operation
            band.set_image_sizes(test_run)
            bands_initsubtract.append(band)

    return bands, bands_initsubtract


def _set_up_directions(parset, bands, field, dry_run=False, test_run=False,
    reset_directions=[]):
    """
    Sets up directions (facets)

    Parameters
    ----------
    parset : dict
        Parset containing processing parameters
    bands : list of Band instances
        Vis data
    field : Direction instance
        Field direction object
    dry_run : bool, optional
        If True, do not run pipelines. All parsets, etc. are made as normal
    test_run : bool, optional
        If True, use test settings. These settings are for testing purposes
        only and will not produce useful results
    reset_directions : list of str, optional
        List of direction names to be reset

    Returns
    -------
    directions : List of Direction instances
        All directions to be used by the run() function
    direction_groups : List of lists of Direction instances
        Groups of directions to be selfcal-ed

    """
    log.info("Building initial sky model...")
    initial_skymodel = factor.directions.make_initial_skymodel(bands[-1])
    log.info('Setting up directions...')

    # First check for user-supplied directions file, then for Factor-generated
    # file from a previous run, then for parameters needed to generate it internally
    dir_parset = parset['direction_specific']
    if 'directions_file' in dir_parset:
        directions = factor.directions.directions_read(dir_parset['directions_file'],
            parset['dir_working'])
    elif os.path.exists(os.path.join(parset['dir_working'], 'factor_directions.txt')):
        directions = factor.directions.directions_read(os.path.join(parset['dir_working'],
            'factor_directions.txt'), parset['dir_working'])
    else:
        if dry_run:
            # Stop here if dry_run is True but no directions file was given
            log.warn('No directions file given. Cannot proceed beyond the '
                'initsubtract operation. Exiting...')
            sys.exit(0)
        elif 'flux_min_jy' not in dir_parset or \
            'size_max_arcmin' not in dir_parset or \
            'separation_max_arcmin' not in dir_parset:
                log.critical('If no directions file is specified, you must '
                    'give values for flux_min_Jy, size_max_arcmin, and '
                    'separation_max_arcmin')
                sys.exit(1)
        else:
            # Make directions from dir-indep sky model of highest-frequency
            # band, as it has the smallest field of view
            log.info("No directions file given. Selecting directions internally...")
            dir_parset['directions_file'] = factor.directions.make_directions_file_from_skymodel(initial_skymodel.copy(),
                dir_parset['flux_min_jy'], dir_parset['size_max_arcmin'],
                dir_parset['separation_max_arcmin'], directions_max_num=dir_parset['max_num'],
                interactive=parset['interactive'])
            directions = factor.directions.directions_read(dir_parset['directions_file'],
                parset['dir_working'])

    # Add the target to the directions list if desired
    target_has_own_facet = dir_parset['target_has_own_facet']
    if target_has_own_facet:
        if 'target_ra' in dir_parset and 'target_dec' in dir_parset and 'target_radius_arcmin' in dir_parset:
            # Make target object
            target = Direction('target', dir_parset['target_ra'], dir_parset['target_dec'],
                factor_working_dir=parset['dir_working'])

            # Check if target is already in directions list. If so, remove it
            nearest, dist = factor.directions.find_nearest(target, directions)
            if dist < dir_parset['target_radius_arcmin']/60.0:
                directions.remove(nearest)

            # Add target to directions list
            directions.append(target)
        else:
            log.critical('target_has_own_facet = True, but target RA, Dec, or radius not found in parset')
            sys.exit(1)

    if 'target_ra' in dir_parset and 'target_dec' in dir_parset and \
        'target_radius_arcmin' in dir_parset:
        target_ra = dir_parset['target_ra']
        target_dec = dir_parset['target_dec']
        target_radius_arcmin = dir_parset['target_radius_arcmin']
    else:
        target_ra = None
        target_dec = None
        target_radius_arcmin = None

    factor.directions.thiessen(directions, s=initial_skymodel,
        check_edges=dir_parset['check_edges'], target_ra=target_ra,
        target_dec=target_dec, target_radius_arcmin=target_radius_arcmin)

    # Warn user if they've specified a direction to reset that does not exist
    direction_names = [d.name for d in directions]
    for name in reset_directions:
        if name not in direction_names and name != 'field':
            log.warn('Direction {} was specified for resetting but does not '
                'exist in current list of directions'.format(name))

    # Make DS9 region files so user can check the facets, etc.
    ds9_facet_reg_file = os.path.join(parset['dir_working'], 'regions', 'facets_ds9.reg')
    factor.directions.make_ds9_region_file(directions, ds9_facet_reg_file)
    ds9_calimage_reg_file = os.path.join(parset['dir_working'], 'regions', 'calimages_ds9.reg')
    factor.directions.make_ds9_calimage_file(directions, ds9_calimage_reg_file)

    # Check with user
    if parset['interactive']:
        print("Facet and DDE calibrator regions saved. Please check that they "
            "are OK before continuing.")
        prompt = "Continue processing (y/n)? "
        answ = raw_input(prompt)
        while answ.lower() not in  ['y', 'n', 'yes', 'no']:
            answ = raw_input(prompt)
        if answ.lower() in ['n', 'no']:
            log.info('Exiting...')
            sys.exit()

    # Select subset of directions to process
    if 'ndir_total' in dir_parset:
        if dir_parset['ndir_total'] > 0 and dir_parset['ndir_total'] <= len(directions):
            directions = directions[:dir_parset['ndir_total']]

    # Set various direction attributes
    log.info("Determining imaging parameters for each direction...")
    for i, direction in enumerate(directions):
        # Set averaging steps and solution intervals for selfcal
        direction.set_averaging_steps_and_solution_intervals(bands[0].chan_width_hz, bands[0].nchan,
            bands[0].timepersample, bands[0].nsamples, len(bands), parset['preaverage_flux_jy'])

        # Set imaging parameters
        direction.set_imaging_parameters(len(bands), parset['wsclean_nbands'],
            initial_skymodel.copy(), test_run=test_run)

        # Set field center
        direction.field_ra = field.ra
        direction.field_dec = field.dec

        # Load previously completed operations (if any)
        direction.load_state()

        # Reset state if specified
        if direction.name in reset_directions:
            direction.do_reset = True
        else:
            direction.do_reset = False

    # Select subset of directions to selfcal
    if target_has_own_facet:
        # Make sure target is not a DDE calibrator
        selfcal_directions = [d for d in directions if d.name != target.name]
    else:
        selfcal_directions = directions
    if 'ndir_selfcal' in dir_parset:
        if dir_parset['ndir_selfcal'] > 0 and dir_parset['ndir_selfcal'] <= len(selfcal_directions):
            selfcal_directions = selfcal_directions[:dir_parset['ndir_selfcal']]

    # Divide directions into groups for selfcal
    direction_groups = factor.directions.group_directions(selfcal_directions,
        n_per_grouping=dir_parset['groupings'], allow_reordering=dir_parset['allow_reordering'])

    # Ensure that target is included in the directions to process if desired
    # (but not for selfcal), just in case it was excluded by one of the cuts
    # above
    if target_has_own_facet:
        names = [d.name for d in directions]
        if target.name not in names:
            directions.append(target)

    return directions, direction_groups
