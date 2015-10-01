"""
Module that preforms the processing
"""
import sys
import os
import shutil
import numpy as np
import logging
import pickle
import lofar.parmdb
from lofarpipe.support.data_map import DataMap
import factor
import factor.directions
import factor.parset
import factor.cluster
from factor.operations.field_ops import *
from factor.operations.facet_ops import *
from factor.lib.scheduler import Scheduler
from factor.lib.direction import Direction


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
    factor._logging.set_level(logging_level)
    log = logging.getLogger('factor')

    parset = factor.parset.parset_read(parset_file)
    parset['logging_level'] = logging_level

    # Set up clusterdesc, node info, scheduler, etc.
    scheduler = _set_up_compute_parameters(parset, log, dry_run)

    # Prepare vis data
    bands, bands_initsubtract = _set_up_bands(parset, log, test_run)

    # Make direction object for the field
    field = Direction('field', bands[0].ra, bands[0].dec,
        factor_working_dir=parset['dir_working'])
    field.set_averaging_steps(bands[0].chan_width_hz, bands[0].nchan,
        bands[0].timepersample)

    # Run initial sky model generation and create empty datasets
    if len(bands_initsubtract) > 0:
        input_bands_full = [b for b in bands_initsubtract if b.skymodel_dirindep is None]
        if len(input_bands_full) > 0:
            log.debug('Running full initial subtract operation for bands: {0}'.
                format([b.name for b in input_bands_full]))
            field = factor.cluster.divide_nodes([field],
                parset['cluster_specific']['node_list'],
                parset['cluster_specific']['ndir_per_node'],
                parset['cluster_specific']['nimg_per_node'],
                parset['cluster_specific']['ncpu'],
                parset['cluster_specific']['fmem'])[0]
            op = InitSubtract(parset, input_bands_full, field)
            scheduler.run(op)
    else:
        log.info("Sky models and SUBTRACTED_DATA_ALL found for all bands. "
            "Skipping initsubtract operation...")

    # Remove bands that failed initsubtract operation
    bands = [b for b in bands if not b.skip]

    # Define directions
    directions, direction_groups = _set_up_directions(parset, bands, field, log,
        dry_run, test_run, reset_directions)

    # Run selfcal and subtract operations on directions
    first_pass = True
    for gindx, direction_group in enumerate(direction_groups):
        log.info('Processing {0} direction(s) in Group {1}'.format(
            len(direction_group), gindx+1))

        # Set up reset of any directions that need it. If the direction has
        # already been through the facetsub operation, we must undo the
        # changes with the facetsubreset operation before we reset
        direction_group_reset = [d for d in direction_group if d.do_reset]
        direction_group_reset_facetsub = [d for d in direction_group_reset if
            'facetsub' in d.completed_operations]
        if len(direction_group_reset_facetsub) > 0:
            direction_group_reset_facetsub = factor.cluster.combine_nodes(
                direction_group_reset_facetsub,
                parset['cluster_specific']['node_list'],
                parset['cluster_specific']['ncpu'],
                parset['cluster_specific']['fmem'])
            ops = [FacetSubReset(parset, d) for d in direction_group_reset_facetsub]
            for op in ops:
                scheduler.run(op)
        for d in direction_group_reset:
            d.reset_state()

        # Divide up the nodes and cores among the directions for the parallel
        # selfcal operations
        direction_group = factor.cluster.divide_nodes(direction_group,
            parset['cluster_specific']['node_list'],
            parset['cluster_specific']['ndir_per_node'],
            parset['cluster_specific']['nimg_per_node'],
            parset['cluster_specific']['ncpu'],
            parset['cluster_specific']['fmem'])

        # Do selfcal on calibrator only
        ops = [FacetSelfcal(parset, bands, d) for d in direction_group]
        scheduler.run(ops)
        if dry_run:
            # For dryrun, skip check
            for d in direction_group:
                d.selfcal_ok = True
        direction_group_ok = [d for d in direction_group if d.selfcal_ok]
        if first_pass:
            if len(direction_group_ok) > 0:
                for d in directions:
                    d.subtracted_data_colname = 'SUBTRACTED_DATA_ALL_NEW'
                first_pass = False

        # Combine the nodes and cores for the serial subtract operations
        direction_group_ok = factor.cluster.combine_nodes(direction_group_ok,
            parset['cluster_specific']['node_list'],
            parset['cluster_specific']['ncpu'],
            parset['cluster_specific']['fmem'])

        # Subtract final model(s) for directions for which selfcal went OK
        ops = [FacetSub(parset, d) for d in direction_group_ok]
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

    # Make final facet images (from final empty datasets) if desired. Also image
    # any facets for which selfcal failed or no selfcal was done
    #
    # TODO: combine new facet sky models and adjust facet edges for new sources
    # (only if all facets are to be re-imaged)
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
        nearest = factor.directions.find_nearest(d, dirs_with_selfcal)
        log.debug('Using solutions from direction {0} for direction {1}.'.format(
            nearest.name, d.name))
        d.dir_dep_parmdb_datamap = nearest.dir_dep_parmdb_datamap
        d.save_state()
    dirs_to_image.extend(dirs_to_transfer)

    if len(dirs_to_image) > 0:
        # Divide up the nodes and cores among the directions for the parallel
        # imaging operations
        dirs_to_image = factor.cluster.divide_nodes(dirs_to_image,
            parset['cluster_specific']['node_list'],
            parset['cluster_specific']['ndir_per_node'],
            parset['cluster_specific']['nimg_per_node'],
            parset['cluster_specific']['ncpu'],
            parset['cluster_specific']['fmem'])

        ops = [FacetImage(parset, bands, d) for d in dirs_to_image]
        scheduler.run(ops)

    # Mosaic the final facet images together
    if parset['make_mosaic']:
        field.facet_image_filenames = []
        field.facet_vertices_filenames = []
        for d in directions:
            facet_image = DataMap.load(d.facet_image_mapfile)[0].file
            field.facet_image_filenames.append(facet_image)
            field.facet_vertices_filenames.append(d.save_file)
        op = MakeMosaic(parset, field)
        scheduler.run(op)

    log.info("Factor has finished :)")


def _set_up_compute_parameters(parset, log, dry_run=False):
    """
    Sets up compute parameters and operation scheduler

    Parameters
    ----------
    parset : dict
        Parset containing processing parameters
    log : logging instance
        Log for output
    dry_run : bool, optional
        If True, do not run pipelines. All parsets, etc. are made as normal

    Returns
    -------
    scheduler : Scheduler instance
        The operation scheduler used by the run() function

    """
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
    ndir_simul = len(parset['cluster_specific']['node_list']) * parset['cluster_specific']['ndir_per_node']
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


def _set_up_bands(parset, log, test_run=False):
    """
    Sets up bands for processing

    Parameters
    ----------
    parset : dict
        Parset containing processing parameters
    log : logging instance
        Log for output
    test_run : bool, optional
        If True, use test settings. These settings are for testing purposes
        only and will not produce useful results

    Returns
    -------
    bands : List of Band instances
        All bands to be used by the run() function
    bands_initsubtract : List of Band instances
        Subset of bands for InitSubtract operation

    """
    bands = []
    from factor.lib.band import Band
    for ms in parset['mss']:
        band = Band(ms, parset['dir_working'], test_run=test_run)

        # Some checks on the dir-indep instrument parmdb
        band.dirindparmdb = os.path.join(band.file, parset['parmdb_name'])
        if parset['parmdb_name'] == 'instrument':
            # Check for special BBS table name
            band.dirindparmdb += '_dirindep'
            if not os.path.exists(band.dirindparmdb):
                if not os.path.exists(os.path.join(band.file, parset['parmdb_name'])):
                    log.critical('Direction-independent instument parmdb not found '
                        'for band {0}'.format(band.file))
                    sys.exit(1)
                log.warn('Direction-independent instument parmdb for band {0} is '
                    'named "instrument". Copying to "instrument_dirindep" so that BBS '
                    'will not overwrite this table...'.format(band.file))
                os.system('cp -r {0} {1}'.format(os.path.join(band.file,
                    parset['parmdb_name']), band.dirindparmdb))
        if not os.path.exists(band.dirindparmdb):
            log.critical('Direction-independent instrument parmdb not found '
                'for band {0}'.format(band.file))
            sys.exit(1)

        try:
            solname = lofar.parmdb.parmdb(band.dirindparmdb).getNames()[0]
        except IndexError:
            log.critical('Direction-independent instument parmdb appears to be empty '
                        'for band {0}'.format(band.file))
            sys.exit(1)

        if 'Real' in solname or 'Imag' in solname:
            # Convert real/imag to phasors
            log.warn('Direction-independent instument parmdb for band {0} contains '
                'real/imaginary values. Converting to phase/amplitude...'.format(band.file))
            band.dirindparmdb = _convert_to_phasors(band.dirindparmdb)
        band.skymodel_dirindep = None
        msbase = os.path.basename(ms)
        if msbase in parset['ms_specific']:
            if 'init_skymodel' in parset['ms_specific'][msbase]:
                band.skymodel_dirindep = parset['ms_specific'][msbase]['init_skymodel']
                if not os.path.exists(band.skymodel_dirindep):
                    log.error('Sky model specified in parset for band {} was '
                        'not found'.format(band.msname))
                    sys.exit(1)
        bands.append(band)

    # Sort bands by frequency
    band_freqs = [band.freq for band in bands]
    bands = np.array(bands)[np.argsort(band_freqs)].tolist()

    # TODO: Check that bands are uniform in number of channels and phase center

    # Determine whether any bands need to be run through the initsubract operation.
    # This operation is needed (only needed if band lacks an initial skymodel or
    # the SUBTRACTED_DATA_ALL column).
    bands_initsubtract = []
    for band in bands:
        if band.skymodel_dirindep is None or not band.has_sub_data:
            # Set image sizes for initsubtract operation
            band.set_image_sizes(test_run)
            bands_initsubtract.append(band)

    return bands, bands_initsubtract


def _convert_to_phasors(real_imag_parmdb_file):
    """
    Converts instrument parmdb from real/imag to phasors

    Parameters
    ----------
    real_imag_parmdb_file : str
        Filename of input parmdb

    Returns
    -------
    phasors_parmdb_file : str
        Filename of ouput phasors parmdb

    """
    phasors_parmdb_file = real_imag_parmdb_file + '_phasors'
    if os.path.exists(phasors_parmdb_file):
        return phasors_parmdb_file

    pdb_in = lofar.parmdb.parmdb(real_imag_parmdb_file)
    pdb_out = lofar.parmdb.parmdb(phasors_parmdb_file, create=True)

    # Get station names
    stations = set([s.split(':')[-1] for s in pdb_in.getNames()])

    # Calculate and store phase and amp values for each station
    parms = pdb_in.getValuesGrid('*')
    for i, s in enumerate(stations):
        if i == 0:
            freqs = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['freqs'])
            freqwidths = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['freqwidths'])
            times = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['times'])
            timewidths = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['timewidths'])

        valIm_00 = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['values'][:, 0])
        valIm_11 = np.copy(parms['Gain:1:1:Imag:{}'.format(s)]['values'][:, 0])
        valRe_00 = np.copy(parms['Gain:0:0:Real:{}'.format(s)]['values'][:, 0])
        valRe_11 = np.copy(parms['Gain:1:1:Real:{}'.format(s)]['values'][:, 0])

        valAmp_00 = np.sqrt((valRe_00**2) + (valIm_00**2))
        valAmp_11 = np.sqrt((valRe_11**2) + (valIm_11**2))
        valPh_00 = np.arctan2(valIm_00, valRe_00)
        valPh_11 = np.arctan2(valIm_11, valRe_11)

        pdb_out.addValues({'Gain:0:0:Phase:{}'.format(s): {'freqs': freqs, 'freqwidths':
            freqwidths, 'times': times, 'timewidths': timewidths, 'values': valPh_00[:,np.newaxis]}})
        pdb_out.addValues({'Gain:1:1:Phase:{}'.format(s): {'freqs': freqs, 'freqwidths':
            freqwidths, 'times': times, 'timewidths': timewidths, 'values': valPh_11[:,np.newaxis]}})
        pdb_out.addValues({'Gain:0:0:Ampl:{}'.format(s): {'freqs': freqs, 'freqwidths':
            freqwidths, 'times': times, 'timewidths': timewidths, 'values': valAmp_00[:,np.newaxis]}})
        pdb_out.addValues({'Gain:1:1:Ampl:{}'.format(s): {'freqs': freqs, 'freqwidths':
            freqwidths, 'times': times, 'timewidths': timewidths, 'values': valAmp_11[:,np.newaxis]}})

    # Write values
    pdb_out.flush()

    return phasors_parmdb_file


def _set_up_directions(parset, bands, field, log, dry_run=False, test_run=False,
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
    log : logging instance
        Log for output
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
    # First check for user-supplied directions file, then for Factor-generated
    # file from a previous run, then for parameters needed to generate it internally
    dir_parset = parset['direction_specific']
    if 'directions_file' in parset:
        directions = factor.directions.directions_read(parset['directions_file'],
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
            # Make directions from dir-indep sky models using flux and size parameters
            log.info("No directions file given. Selecting directions internally...")
            parset['directions_file'] = factor.directions.make_directions_file_from_skymodel(bands,
                dir_parset['flux_min_jy'], dir_parset['size_max_arcmin'],
                dir_parset['separation_max_arcmin'], directions_max_num=dir_parset['max_num'],
                interactive=parset['interactive'])
            directions = factor.directions.directions_read(parset['directions_file'],
                parset['dir_working'])

    # Add the target to the directions list if desired
    target_has_own_facet = dir_parset['target_has_own_facet']
    if target_has_own_facet:
        if 'target_ra' in dir_parset and 'target_dec' in dir_parset and 'target_radius_arcmin' in dir_parset:
            # Make target object
            target = Direction('target', dir_parset['target_ra'], dir_parset['target_dec'],
                factor_working_dir=parset['dir_working'])

            # Check if target is already in directions list. If so, remove it
            nearest = factor.directions.find_nearest(target, directions)
            dist = factor.directions.calculateSeparation(target.ra, target.dec,
                nearest.ra, nearest.dec)
            if dist.value < dir_parset['target_radius_arcmin']/60.0:
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

    factor.directions.thiessen(directions, band=bands[0],
        check_edges=dir_parset['check_edges'], target_ra=target_ra,
        target_dec=target_dec, target_radius_arcmin=target_radius_arcmin)

    # Set various direction attributes
    for i, direction in enumerate(directions):
        # Set averaging steps
        direction.set_averaging_steps(bands[0].chan_width_hz, bands[0].nchan,
            bands[0].timepersample)

        # Set image sizes
        direction.set_image_sizes(test_run=test_run)

        # Set number of bands and channels for images
        direction.nbands = len(bands)
        if direction.nbands > 5:
            direction.nchannels = int(round(float(direction.nbands)/
                float(parset['wsclean_nbands'])))
        else:
            direction.nchannels = 1

        # Set field center
        direction.field_ra = field.ra
        direction.field_dec = field.dec

        # Load previous state (if any)
        direction.load_state()

        # Reset state if specified
        if direction.name in reset_directions:
            direction.do_reset = True
        else:
            direction.do_reset = False

    # Warn user if they've specified a direction to reset that does not exist
    direction_names = [d.name for d in directions]
    for name in reset_directions:
        if name not in direction_names:
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
    if 'ndir_total' in parset['direction_specific']:
        if parset['direction_specific']['ndir_total'] > 0 and \
            parset['direction_specific']['ndir_total'] <= len(directions):
            directions = directions[:parset['direction_specific']['ndir_total']]

    # Select subset of directions to selfcal
    if target_has_own_facet:
        # Make sure target is not a DDE calibrator
        selfcal_directions = [d for d in directions if d.name != target.name]
    else:
        selfcal_directions = directions
    if 'ndir_selfcal' in parset['direction_specific']:
        if parset['direction_specific']['ndir_selfcal'] > 0 and \
            parset['direction_specific']['ndir_selfcal'] <= len(selfcal_directions):
            selfcal_directions = selfcal_directions[:parset['direction_specific']['ndir_selfcal']]

    # Divide directions into groups for selfcal
    direction_groups = factor.directions.group_directions(selfcal_directions,
        one_at_a_time=parset['direction_specific']['one_at_a_time'],
        n_per_grouping=parset['direction_specific']['groupings'],
        allow_reordering=parset['direction_specific']['allow_reordering'])

    # Ensure that target is included in the directions to process if desired
    # (but not for selfcal), just in case it was excluded by one of the cuts
    # above
    if target_has_own_facet:
        names = [d.name for d in directions]
        if target.name not in names:
            directions.append(target)

    return directions, direction_groups
