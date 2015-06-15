"""
Module that holds the process run() function called by runfactor
"""
import sys
import os
import numpy as np
import logging
import pickle
import factor
import factor.directions
import factor.parset
import factor.cluster
from factor.operations.field_ops import *
from factor.operations.facet_ops import *
from factor.lib.scheduler_mp import Scheduler
from factor.lib.direction import Direction


def run(parset_file, logging_level='info', dry_run=False, test_run=False):
    """
    Processes a dataset

    Parameters
    ----------
    parset_file : str
        Filename of parset containing processing parameters
    logging_level : str, optional
        One of 'degug', 'info', 'warning'
    dry_run : bool, optional
        If True, do not run pipelines
    test_run : bool, optional
        If True, use test settings
    """
    factor._logging.set_level(logging_level)
    log = logging.getLogger('factor')

    parset = factor.parset.parset_read(parset_file)
    parset['logging_level'] = logging_level

    # Prepare vis data
    bands = []
    from factor.lib.band import Band
    for ms in parset['mss']:
        band = Band(ms, parset['dir_working'])
        band.dirindparmdb = os.path.join(band.file, parset['parmdb_name'])
        if not os.path.exists(band.dirindparmdb):
            log.critical('Direction-independent instument parmdb not found '
                'for band {0}'.format(band.file))
            sys.exit(1)
        band.skymodel_dirindep = None
        msbase = os.path.basename(ms)
        if msbase in parset['ms_specific']:
            if 'init_skymodel' in parset['ms_specific'][msbase]:
                band.skymodel_dirindep = parset['ms_specific'][msbase]['init_skymodel']
        bands.append(band)

    # Get clusterdesc, node info, etc.
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
    factor.cluster.find_executables(parset)

    # Set up scheduler for operations (pipeline runs)
    scheduler = Scheduler(max_procs=len(parset['cluster_specific']['node_list']),
        dry_run=dry_run)

    # Make direction object for the field
    field = Direction('field', bands[0].ra, bands[0].dec,
        factor_working_dir=parset['dir_working'])
    field.imsize_high_res = 6144 # TODO: calculate image size for each field
    field.imsize_low_res = 4800 # TODO: calculate image size for each field

    # Run initial sky model generation and create empty datasets. First check that
    # this operation is needed (only needed if band lacks an initial skymodel or
    # the SUBTRACTED_DATA_ALL column).
    bands_init_subtract = []
    for band in bands:
        if band.skymodel_dirindep is None or not band.has_sub_data:
            bands_init_subtract.append(band)
    if len(bands_init_subtract) > 0:
        op = InitSubtract(parset, bands_init_subtract, field)
        scheduler.run(op)
        field.cleanup()
    else:
        log.info("Sky models found for all MS files. Skipping initial subtraction "
            "operation")

    # Define directions. First check for user-supplied file, then for Factor-generated
    # file from a previous run, then for parameters needed to generate it internally
    if 'directions_file' in parset:
        directions = factor.directions.directions_read(parset['directions_file'],
            parset['dir_working'])
    elif os.path.exists(os.path.join(parset['dir_working'], 'factor_directions.txt')):
        directions = factor.directions.directions_read(os.path.join(parset['dir_working'],
            'factor_directions.txt'), parset['dir_working'])
    else:
        dir_parset = parset['direction_specific']
        if 'flux_min_jy' not in dir_parset or \
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
                parset['direction_specific']['flux_min_jy'],
                parset['direction_specific']['size_max_arcmin'],
                parset['direction_specific']['separation_max_arcmin'],
                directions_max_num=parset['direction_specific']['max_num'],
                interactive=parset['interactive'])
            directions = factor.directions.directions_read(parset['directions_file'],
                parset['dir_working'])

    # Load polygons from previous run if possible
    polys_file = os.path.join(parset['dir_working'], 'regions', 'factor_facets.pkl')
    if os.path.exists(polys_file):
        with open(polys_file, 'r') as f:
            polys, widths = pickle.load(f)
            widths = [w[0] for w in widths]
    else:
        polys, widths = factor.directions.thiessen(directions,
            check_edges=parset['direction_specific']['check_edges'])
        with open(polys_file, 'wb') as f:
            pickle.dump([polys, widths], f)

    # Set various direction attributes
    for i, direction in enumerate(directions):
        exists = direction.load_state()
        if not exists:
            direction.vertices = polys[i]
            direction.width = widths[i]

            # Set image sizes
            direction.facet_imsize = getOptimumSize(direction.width * 3600.0 / 1.5
                * 1.15) # full facet has 15% padding to avoid aliasing issues with ft
            direction.cal_imsize = getOptimumSize(direction.cal_radius_deg * 3600.0
                / 1.5 * 1.5) # cal size has 50% padding

            # Make CASA region files for use during clean
            reg_file = os.path.join(parset['dir_working'], 'regions', direction.name+'.rgn')
            factor.directions.make_region_file(direction.vertices, reg_file)
            direction.reg = reg_file

            # Set number of bands and channels
            direction.nbands = len(bands)
            direction.nchannels = np.int(np.ceil(np.float(direction.nbands/np.float(5))))

            # Set field center
            direction.field_ra = field.ra
            direction.field_dec = field.dec

            # Save direction state
            direction.save_state()

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
            self.log.info('Exiting...')
            sys.exit()

    # Select subset of directions to process
    if 'ndir' in parset['direction_specific']:
        if parset['direction_specific']['ndir'] > 0 and \
            parset['direction_specific']['ndir'] <= len(directions):
            directions = directions[:parset['direction_specific']['ndir']]

    direction_groups = factor.directions.group_directions(directions,
        one_at_a_time=parset['direction_specific']['one_at_a_time'],
        n_per_grouping=parset['direction_specific']['groupings'])

    # Iterate over direction groups
    first_pass = True
    for direction_group in direction_groups:
        log.info('Processing {0} direction(s) in parallel in this group'.format(
            len(direction_group)))

        # Divide up the nodes among the directions
        node_list = parset['cluster_specific']['node_list']
        if len(direction_group) >= len(node_list):
            for i in range(len(direction_group)-len(node_list)):
                node_list.append(node_list[i])
            hosts = [[n] for n in node_list]
        else:
            parts = len(direction_group)
            hosts = [node_list[i*len(node_list)//parts:
                (i+1)*len(node_list)//parts] for i in range(parts)]
        for d, h in zip(direction_group, hosts):
            d.hosts = h
            d.save_state()

        # Add calibrator(s) to empty datasets. These operations
        # must be done in series
        ops = [FacetAdd(parset, bands, d) for d in direction_group]
        for op in ops:
            scheduler.run(op)

        # Do selfcal on calibrator only
        ops = [FacetSelfcal(parset, bands, d) for d in direction_group]
        scheduler.run(ops)

        # Subtract final model(s) from empty field datasets. These operations
        # must be done in series and only on the directions that passed the
        # selfcal check. Also, after this operation is complete for any
        # direction, set flag to indicate all subsequent directions should use
        # new subtracted data column
        if dry_run:
            # For dryrun, skip check
            for d in direction_group:
                d.selfcal_ok = True
                d.save_state()
        direction_group_ok = [d for d in direction_group if d.selfcal_ok]
        if first_pass:
            if len(direction_group_ok) > 0:
                # Only use new data if at least one direction is OK
                for i, d in enumerate(directions):
                    # Set flag for *all* directions except first one
                    if i > 0:
                        d.use_new_sub_data = True
                first_pass = False
        else:
            for d in direction_group_ok:
                d.use_new_sub_data = True
        ops = [FacetSub(parset, bands, d) for d in direction_group if d.selfcal_ok]
        for op in ops:
            scheduler.run(op)

        # Lastly, stop Factor if selfcal for any direction in this group failed
        for d in direction_group:
            all_good = True
            if not d.selfcal_ok:
                log.error('Selfcal failed for direction {0}. Please check '
                    'the settings for this direction.'.format(d.name))
                if parset['interactive']:
                    prompt = "Continue with this direction anyway (y/n)? "
                    answ = raw_input(prompt)
                    while answ.lower() not in  ['y', 'n', 'yes', 'no']:
                        answ = raw_input(prompt)
                    if answ.lower() in ['n', 'no']:
                        d.reset_state()
                        all_good = False
                    else:
                        d.selfcal_ok = True
                        d.save_state()
                else:
                    d.reset_state()
                    all_good = False
            else:
                d.save_state()
        if not all_good:
            self.log.info('Exiting...')
            sys.exit(1)

        # Clean up files
        for d in direction_group:
            d.cleanup()

    # Make final facet images (from final empty datasets) if desired
    dirs_to_image = [d for d in directions if d.make_final_image]
    if len(dirs_to_image) > 0:
        ops = [FacetAddAllFinal(parset, bands, d) for d in dirs_to_image]
        for op in ops:
            scheduler.run(op)
        ops = [FacetImageFinal(parset, bands, d) for d in dirs_to_image]
        scheduler.run(ops)

    # Mosaic the final facet images together
    if parset['make_mosaic']:
        op = MakeMosaic(parset, directions)
        scheduler.run(op)

    log.info("Factor has finished :)")


def getOptimumSize(size):
    """
    Gets the nearest optimum image size

    Taken from the casa source code (cleanhelper.py)

    Parameters
    ----------
    size : int
        Target image size in pixels

    Returns
    -------
    optimum_size : int
        Optimum image size nearest to target size

    """
    import numpy

    def prime_factors(n, douniq=True):
        """ Return the prime factors of the given number. """
        factors = []
        lastresult = n
        sqlast=int(numpy.sqrt(n))+1
        if n == 1:
            return [1]
        c=2
        while 1:
             if (lastresult == 1) or (c > sqlast):
                 break
             sqlast=int(numpy.sqrt(lastresult))+1
             while 1:
                 if(c > sqlast):
                     c=lastresult
                     break
                 if lastresult % c == 0:
                     break
                 c += 1

             factors.append(c)
             lastresult /= c

        if (factors==[]): factors=[n]
        return  numpy.unique(factors).tolist() if douniq else factors

    n = int(size)
    if (n%2 != 0):
        n+=1
    fac=prime_factors(n, False)
    for k in range(len(fac)):
        if (fac[k] > 7):
            val=fac[k]
            while (numpy.max(prime_factors(val)) > 7):
                val +=1
            fac[k]=val
    newlarge=numpy.product(fac)
    for k in range(n, newlarge, 2):
        if ((numpy.max(prime_factors(k)) < 8)):
            return k
    return newlarge
