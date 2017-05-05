"""
Module that holds all archiving functions
"""
import os
import shutil
import logging
import sys
import numpy as np
import subprocess
import factor._logging
import factor.parset
import factor.directions
from factor.lib.direction import Direction
from lofarpipe.support.data_map import DataMap
from factor.scripts import sort_times_into_freqGroups
from lofarpipe.support.utilities import create_directory
import glob
import pickle

log = logging.getLogger('factor:archive')


def check_existing_files(mapfile):
    """
    Checks if files in input mapfile exist

    Parameters
    ----------
    mapfile : str
        Filename of mapfile to check

    Returns
    -------
    file : list
        List of files

    """

    all_exist = True
    all_files = []
    log.info('Checking for existing files...')
    try:
        datamap = DataMap.load(mapfile)
        for item in datamap:
            # Handle case in which item.file is a Python list
            if item.file[0] == '[' and item.file[-1] == ']':
                files = item.file.strip('[]').split(',')
            else:
                files = [item.file]
            for f in files:
                if not os.path.exists(f):
                    all_exist = False
            all_files.extend(files)
        if all_exist:
            log.info('...all files exist')
        else:
            log.warning('...one or more files not found')
        return all_files
    except IOError:
        return []


def load_directions(parset_file):
    """
    Return directions for a run
    """
    # Read parset
    orig_dir = os.path.abspath('.')
    parset = factor.parset.parset_read(parset_file, use_log_file=False)
    os.chdir(orig_dir)

    # Load directions. First check for user-supplied directions file then for
    # Factor-generated file from a previous run
    direction_list = []
    dir_parset = parset['direction_specific']
    if 'directions_file' in dir_parset:
        directions = factor.directions.directions_read(dir_parset['directions_file'],
            parset['dir_working'])
    elif os.path.exists(os.path.join(parset['dir_working'], 'factor_directions.txt')):
        directions = factor.directions.directions_read(os.path.join(parset['dir_working'],
            'factor_directions.txt'), parset['dir_working'])
    else:
        log.error('No directions found. Please run this tool after '
            'the directions have been defined')
        sys.exit(1)

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

            # Add target to directions list
            directions.append(target)
        else:
            log.critical('target_has_own_facet = True, but target RA, Dec, or radius not found in parset')
            sys.exit(1)

    for direction in directions:
        has_state = direction.load_state()
        if has_state:
            direction_list.append(direction)

    return direction_list, parset


def copy(path_from, dir_to, clobber, use_symlinks=False):
    """
    Copy a file or directory

    Parameters
    ----------
    path_from : str
        Input file or directory
    dir_to : str
        Output directory
    clobber : bool
        Clobber existing file or directory?
    use_symlinks : bool, optional
        Use symlinks instead of copying files?

    """
    if not os.path.exists(path_from):
        log.warning('{} not found. Please check the '
            'working directory'.format(path_from))
        return

    path_to = os.path.join(dir_to, os.path.basename(path_from))
    if os.path.exists(path_to):
        if not clobber:
            log.warning(' Destination "{}" exists and clobber = False. '
                'Skipping it...'.format(path_to))
            return
    else:
        create_directory(dir_to)

    if use_symlinks:
        if os.path.exists(path_to):
            p = subprocess.Popen('rm -rf {0}'.format(path_to), shell=True,
                stdout=subprocess.PIPE)
            r = p.communicate()
        os.symlink(path_from, path_to)
    else:
        p = subprocess.Popen('rsync -a {0} {1}'.format(path_from, dir_to),
            shell=True, stdout=subprocess.PIPE)
        r = p.communicate()
        if p.returncode != 0:
            log.critical('rsync exited abnormally when attempting to archive {}'.format(path_from))
            sys.exit(1)


def dppp_concat(mslist, msout):
    """
    Run DPPP to concat a list of files

    Parameters
    ----------
    mslist : str
        List of input ms files, given as a string (e.g., '[ms1,ms2,...]')
    msout : str
        Filename of output ms file

    """
    # Call DPPP
    p = subprocess.Popen("DPPP msin={0} steps=[] msout={1} msin.missingdata=true "
        "msin.orderms=false".format(mslist, msout), shell=True, stdout=subprocess.PIPE)
    r = p.communicate()

    if p.returncode != 0:
        log.critical('DPPP exited abnormally when attempting to concat {}'.format(mslist))
        sys.exit(1)


def archive(parset_file, directions, dir_output, full=False, archive_subdata=False,
    archive_state=False, archive_misc=True, archive_images=True,
    archive_inst=False, archive_pipestate=False, archive_models=False,
    archive_plots=True, clobber=False):
    """
    Archives data from a Factor run

    Parameters
    ----------
    parset_file : str
        Filename of Factor parset for run of interest
    directions : list of str
        List of direction names for which to archive the calibrated data
    dir_output : str
        Name of output directory where archived data will be stored
    full : bool, optional
        Make a full archive suitable for resuming?
    archive_subdata : bool, optional
        Archive the subtracted data MS files?
    archive_state : bool, optional
        Archive the state files?
    archive_misc : bool, optional
        Archive miscelaneous files?
    archive_images : bool, optional
        Archive the facet and field images?
    archive_inst : bool, optional
        Archive the instrument tables?
    archive_pipestate : bool, optional
        Archive the pipeline state files?
    archive_models : bool, optional
        Archive the sky models?
    archive_plots : bool, optional
        Archive the selfcal plots?
    clobber : bool, optional
        Clobber existing files in output directory?

    """
    # Read in parset and get directions
    all_directions, parset = load_directions(parset_file)
    if len(all_directions) == 0:
        log.error('No directions found in Factor working directory. Please check '
            'the parset')
        sys.exit(1)
    all_names = [d.name for d in all_directions]
    if len(directions) != 0:
        if directions[0].lower() == 'all':
            directions = all_names
        for dname in directions:
            if dname not in all_names:
                log.warning('Direction {} not found. Skipping it...'.format(dname))

    if full:
        # Archive everything
        archive_subdata = True
        archive_state = True
        archive_misc = True
        archive_images = True
        archive_inst = True
        archive_pipestate = True
        archive_models = True
        archive_plots = True

    working_dir = all_directions[0].working_dir
    if archive_subdata:
        log.info('Archiving subtracted data files...')
        chunks_dir = os.path.join(working_dir, 'chunks')
        copy(chunks_dir, dir_output, clobber)

    if archive_state:
        log.info('Archiving state files...')
        state_dir = os.path.join(working_dir, 'state')
        copy(state_dir, dir_output, clobber)

    if archive_misc:
        log.info('Archiving miscelaneous files...')
        misc_dir = os.path.join(dir_output, 'misc')
        if 'directions_file' in parset['direction_specific']:
            directions_file = parset['direction_specific']['directions_file']
        else:
            directions_file = os.path.join(working_dir, 'factor_directions.txt')
        file_list = [directions_file,
                     parset_file,
                     '{}/factor.log'.format(working_dir),
                     '{}/regions/facets_ds9.reg'.format(working_dir),
                     '{}/regions/calimages_ds9.reg'.format(working_dir)]
        for f in file_list:
            copy(f, misc_dir, clobber)

    if archive_images:
        log.info('Archiving field images...')
        file_list = glob.glob(os.path.join(working_dir, 'results',
            'field*', 'field', '*.fits'))
        if len(file_list) == 0:
            log.warning('No field images found.')
        else:
            for i, f in enumerate(file_list):
                log.info('  Archiving image {0} of {1}...'.format(i+1, len(file_list)))
                subdir = f.split('/')[-3]
                image_dir = os.path.join(dir_output, 'images', 'field', subdir)
                copy(f, image_dir, clobber)

    if archive_models:
        log.info('Archiving direction-independent sky models...')
        band_state_files = glob.glob(os.path.join(working_dir, 'state',
            'Band_*'))
        file_list = []
        band_list = []
        for bf in band_state_files:
            try:
                with open(bf, 'r') as f:
                    b = pickle.load(f)
                    file_list.append(b['skymodel_dirindep'])
                    band_list.append(b['name'])
            except:
                pass
        for i, f in enumerate(file_list):
            skymodel_dir = os.path.join(dir_output, 'chunks', band_list[i])
            log.info('  Copying sky model file {0} of {1}...'.format(i+1, len(file_list)))
            copy(f, skymodel_dir, clobber)

    for d in all_directions:
        if archive_images:
            log.info('Archiving facet images for direction {}...'.format(d.name))
            file_list = glob.glob(os.path.join(working_dir, 'results',
                'facetimage*', d.name, '*full2*image.fits'))
            if len(file_list) == 0:
                log.warning('No facet images found for direction {}.'.format(d.name))
            else:
                for i, f in enumerate(file_list):
                    subdir = f.split('/')[-3]
                    image_dir = os.path.join(dir_output, 'images', d.name, subdir)
                    copy(f, image_dir, clobber)

        if archive_models:
            log.info('Archiving sky models for direction {}...'.format(d.name))
            if hasattr(d, 'sourcedb_new_facet_sources'):
                file_list = check_existing_files(d.sourcedb_new_facet_sources)
            else:
                file_list = []
            if len(file_list) == 0:
                log.warning('No sky models found for direction {}.'.format(d.name))
            else:
                sourcedb_dir = os.path.join(dir_output, 'sky_models', d.name)
                for i, f in enumerate(file_list):
                    log.info('  Copying sky model file {0} of {1}...'.format(i+1, len(file_list)))
                    copy(f, sourcedb_dir, clobber)

        if archive_inst:
            log.info('Archiving instrument tables for direction {}...'.format(d.name))
            if hasattr(d, 'converted_parmdb_mapfile'):
                file_list = check_existing_files(d.converted_parmdb_mapfile)
            else:
                file_list = []
            if hasattr(d, 'preapply_parmdb_mapfile'):
                file_list.append(check_existing_files(d.preapply_parmdb_mapfile))
            if len(file_list) == 0:
                log.warning('No instrument tables found for direction {}.'.format(d.name))
            else:
                inst_table_dir = os.path.join(dir_output, 'instrument_tables', d.name)
                for i, f in enumerate(file_list):
                    log.info('  Copying instrument table file {0} of {1}...'.format(i+1, len(file_list)))
                    copy(f, inst_table_dir, clobber)

        if archive_plots:
            log.info('Archiving plots for direction {}...'.format(d.name))
            file_list = glob.glob(os.path.join(working_dir, 'results', 'facetselfcal', d.name, '*png'))
            if len(file_list) == 0:
                file_list = glob.glob(os.path.join(working_dir, 'results', 'facetpeel', d.name, '*png'))
            if len(file_list) == 0:
                file_list = glob.glob(os.path.join(working_dir, 'results', 'outlierpeel', d.name, '*png'))
            if len(file_list) == 0:
                log.warning('No plots found for direction {}.'.format(d.name))
            else:
                plot_dir = os.path.join(dir_output, 'plots', d.name)
                for i, f in enumerate(file_list):
                    copy(f, plot_dir, clobber)

        if archive_pipestate:
            log.info('Archiving pipeline state files for direction {}...'.format(d.name))
            file_list = glob.glob(os.path.join(working_dir, 'results', 'facetselfcal', d.name, 'mapfiles', '*'))
            op_name = 'facetselfcal'
            if len(file_list) == 0:
                file_list = glob.glob(os.path.join(working_dir, 'results', 'facetpeel', d.name, 'mapfiles', '*'))
                op_name = 'facetpeel'
            if len(file_list) == 0:
                file_list = glob.glob(os.path.join(working_dir, 'results', 'outlierpeel', d.name, 'mapfiles', '*'))
                op_name = 'outlierpeel'
            if len(file_list) == 0:
                log.warning('No pipeline state files found for direction {}.'.format(d.name))
            else:
                mapfile_dir = os.path.join(dir_output, 'pipeline_state', d.name, op_name)
                for f in file_list:
                    copy(f, mapfile_dir, clobber)

            # Also archive "final_image" mapfile for facetimage (needed for mosaicking)
            file_list = glob.glob(os.path.join(working_dir, 'results',
                'facetimage*', d.name, 'mapfiles', 'final_image.mapfile'))
            if len(file_list) > 0:
                for i, f in enumerate(file_list):
                    subdir = f.split('/')[-4]
                    mapfile_dir = os.path.join(dir_output, 'pipeline_state', d.name, subdir)
                    copy(f, mapfile_dir, clobber)

        if d.name in directions:
            log.info('Archiving calibrated data for direction {}...'.format(d.name))
            if hasattr(d, 'image_data_mapfile'):
                file_list = check_existing_files(d.image_data_mapfile)
            else:
                file_list = []
            if len(file_list) == 0:
                log.warning('No data found for direction {}. Skipping it...'.format(d.name))
                continue

            # Make the output directory
            cal_data_dir = os.path.join(dir_output, 'calibrated_data', d.name)
            create_directory(cal_data_dir)

            # Sort the files into time chunks
            data_mapfile = d.name+'_calibrated_data.mapfile'
            sort_times_into_freqGroups.main(file_list, filename=data_mapfile,
                mapfile_dir=cal_data_dir)

            # Read the new, grouped file lists
            datamap = DataMap.load(os.path.join(cal_data_dir, data_mapfile))

            # Run DPPP to concatenate each time chunk in frequency
            nchunks = len(datamap)
            for i, item in enumerate(datamap):
                log.info('  Concatenating files for time chunk {0} of {1}...'.format(i+1, nchunks))
                outfile = os.path.join(cal_data_dir, '{0}_calibrated_data_chunk{1}.ms'.format(d.name, i))
                if os.path.exists(outfile):
                    if not clobber:
                        log.warning(' Output file for this chuck exists and clobber = False. Skipping it...')
                        continue
                    else:
                        os.system('rm -rf {0}'.format(outfile))
                dppp_concat(item.file, outfile)

            # Clean up
            os.system('rm -f {0}'.format(os.path.join(cal_data_dir, data_mapfile)))
            os.system('rm -f {0}_groups'.format(os.path.join(cal_data_dir, data_mapfile)))

    log.info('Archiving complete.')
