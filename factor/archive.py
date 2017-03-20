"""
Module that holds all archiving functions
"""
import os
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

    return direction_list


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


def archive(parset_file, directions, dir_output, archive_subdata=False, clobber=False):
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
    archive_subdata : bool, optional
        Archive the subtracted data MS files?
    clobber : bool, optional
        Clobber existing files in output directory?

    """
    # Read in parset and get directions
    all_directions = load_directions(parset_file)
    if len(all_directions) == 0:
        log.error('No directions found in Factor working directory. Please check '
            'the parset')
        sys.exit(1)
    all_names = [d.name for d in all_directions]
    if directions[0].lower() == 'all':
        directions = all_names
    for dname in directions:
        if dname not in all_names:
            log.warning('Direction {} not found. Skipping it...'.format(dname))

    if archive_subdata:
        # Copy chunks directory
        chunks_dir = os.path.join(all_directions[0].working_dir, 'chunks')
        if not os.path.exists(chunks_dir):
            log.critical('No subtracted data files found. Please check the '
                'working directory')
            sys.exit(1)
        os.system('cp -r {0} {1}'.format(chunks_dir, os.path.join(dir_output, '.')))

    # Copy the field images
    log.info('Archiving field images...')
    file_list = glob.glob(os.path.join(all_directions[0].working_dir, 'results',
        'field*', 'field', '*.fits'))
    if len(file_list) == 0:
        log.warning('No field images found.')
    else:
        for i, f in enumerate(file_list):
            log.info('  Copying image {0} of {1}...'.format(i+1, len(file_list)))
            dirs = f.split('/')
            for d in dirs:
                if 'fieldmosaic' in d:
                    subdir = d
                    break
            image_dir = os.path.join(dir_output, 'images', subdir)
            create_directory(image_dir)
            outfile = os.path.join(image_dir, os.path.basename(f))
            if os.path.exists(outfile):
                if not clobber:
                    log.warning(' Output file for this image exists and clobber = False. Skipping it...')
                    continue
                else:
                    os.system('rm -rf {0}'.format(outfile))
            os.system('cp -r {0} {1}'.format(f, outfile))

    for d in all_directions:
        log.info('Archiving sky models for direction {}...'.format(d.name))
        if hasattr(d, 'sourcedb_new_facet_sources'):
            file_list = check_existing_files(d.sourcedb_new_facet_sources)
        else:
            file_list = []
        if len(file_list) == 0:
            log.warning('No sky models found for direction {}.'.format(d.name))
        else:
            sourcedb_dir = os.path.join(dir_output, 'sky_models', d.name)
            create_directory(sourcedb_dir)
            for i, f in enumerate(file_list):
                log.info('  Copying sky model file {0} of {1}...'.format(i+1, len(file_list)))
                outfile = os.path.join(sourcedb_dir, os.path.basename(f))
                if os.path.exists(outfile):
                    if not clobber:
                        log.warning(' Output file for this sky model exists and clobber = False. Skipping it...')
                        continue
                    else:
                        os.system('rm -rf {0}'.format(outfile))
                os.system('cp -r {0} {1}'.format(f, outfile))

        log.info('Archiving instrument tables for direction {}...'.format(d.name))
        if hasattr(d, 'converted_parmdb_mapfile'):
            file_list = check_existing_files(d.converted_parmdb_mapfile)
        else:
            file_list = []
        if len(file_list) == 0:
            log.warning('No instrument tables found for direction {}.'.format(d.name))
        else:
            inst_table_dir = os.path.join(dir_output, 'instrument_tables', d.name)
            create_directory(inst_table_dir)
            for i, f in enumerate(file_list):
                log.info('  Copying instrument table file {0} of {1}...'.format(i+1, len(file_list)))
                outfile = os.path.join(inst_table_dir, os.path.basename(f))
                if os.path.exists(outfile):
                    if not clobber:
                        log.warning(' Output file for this instrument table exists and clobber = False. Skipping it...')
                        continue
                    else:
                        os.system('rm -rf {0}'.format(outfile))
                os.system('cp -r {0} {1}'.format(f, outfile))

        log.info('Archiving plots for direction {}...'.format(d.name))
        file_list = glob.glob(os.path.join(d.working_dir, 'results', 'facetselfcal', d.name, '*png'))
        if len(file_list) == 0:
            log.warning('No plots found for direction {}.'.format(d.name))
        else:
            plot_dir = os.path.join(dir_output, 'plots', d.name)
            create_directory(plot_dir)
            for i, f in enumerate(file_list):
                outfile = os.path.join(plot_dir, os.path.basename(f))
                os.system('cp -r {0} {1}'.format(f, outfile))

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
