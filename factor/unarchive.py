"""
Module that holds all unarchiving functions
"""
import os
import shutil
import logging
import sys
import numpy as np
import subprocess
import factor._logging
from factor.archive import copy
from lofarpipe.support.data_map import DataMap, DataProduct
import glob
import pickle

log = logging.getLogger('factor:unarchive')


def update_state(dir_input):
    """
    Updates the paths in mapfiles or state files

    Parameters
    ----------
    dir_input : str
        Directory containing files to update

    """
    file_list = glob.glob(os.path.join(dir_input, '*'))

    if dir_input.endswith('mapfiles'):
        # Assume path is a pipeline mapfiles directory. In this case, we can
        # simply substitute the new working_dir for the old one in each of the
        # mapfiles
        working_dir = dir_input.split('/results/')[0]
        for f in file_list:
            map = DataMap.load(f)
            for item in map:
                if '/' in item.file:
                    old_working_dir = item.file.split('/results/')[0]
                    item.file = item.file.replace(old_working_dir, working_dir)
            map.save(f)
    elif dir_input.endswith('state'):
        # Assume path is the Factor state directory. In this case, we can try to
        # load files as pickled state files and look for paths inside. If found,
        # substitute new working_dir for the old one
        working_dir = os.path.dirname(dir_input)
        for f in file_list:
            try:
                with open(f, "rb") as fp:
                    d = pickle.load(fp)
                    for k, v in d.iteritems():
                        if type(v) is str:
                            if k == 'working_dir':
                                d[k] = working_dir
                            if '/' in v:
                                for infix in ['results/', 'state/', 'chunks/']:
                                    parts = v.split(infix)
                                    if len(parts) > 1:
                                        d[k] = os.path.join(working_dir, infix, parts[-1])
                        elif type(v) is list:
                            for i, l in enumerate(v):
                                if '/' in l:
                                    for infix in ['results/', 'state/', 'chunks/']:
                                        parts = l.split(infix)
                                        if len(parts) > 1:
                                            v[i] = os.path.join(working_dir, infix, parts[-1])
                            d[k] = v
                with open(f, "w") as fp:
                    pickle.dump(d, fp)
            except:
                pass


def unarchive(dir_input, dir_output, use_symlinks=False, clobber=False):
    """
    Unarchives data from a Factor archive

    Parameters
    ----------
    dir_input : str
        Path to Factor archive
    dir_output : str
        Name of output directory
    use_symlinks : bool, optional
        Use symlinks for subtracted data files instead of copying them?
    clobber : bool, optional
        Clobber existing files in output directory?

    """
    dir_input = os.path.abspath(dir_input)
    dir_output = os.path.abspath(dir_output)

    log.info('Unarchiving subtracted data files...')
    chunks_dir = os.path.join(dir_input, 'chunks')
    copy(chunks_dir, dir_output, clobber, use_symlinks=use_symlinks)

    log.info('Unarchiving state files...')
    state_dir = os.path.join(dir_input, 'state')
    copy(state_dir, dir_output, clobber)
    update_state(os.path.join(dir_output, 'state'))

    # Derive direction names from directories
    direction_dirs = glob.glob(os.path.join(dir_input, 'sky_models', '*'))
    direction_names = [os.path.basename(d) for d in direction_dirs]

    for d in direction_names:
        log.info('Unarchiving facet images for direction {}...'.format(d))
        direction_dir = os.path.join(dir_input, 'images', d)
        file_list = glob.glob(os.path.join(dir_input, 'results',
            'facetimage*', d, '*full2*image.fits'))
        if len(file_list) == 0:
            log.warning('No facet images found for direction {}.'.format(d))
        else:
            for i, f in enumerate(file_list):
                subdir = f.split('/')[-3]
                image_dir = os.path.join(dir_output, 'results', subdir, d)
                copy(f, image_dir, clobber)

        log.info('Unarchiving sky models for direction {}...'.format(d))
        direction_dir = os.path.join(dir_input, 'sky_models', d)
        file_list = glob.glob(os.path.join(direction_dir, '*'))
        sourcedb_dir = os.path.join(dir_output, 'results', 'facetselfcal', d)
        for i, f in enumerate(file_list):
            log.info('  Copying sky model file {0} of {1}...'.format(i+1, len(file_list)))
            copy(f, sourcedb_dir, clobber)

        log.info('Unarchiving instrument tables for direction {}...'.format(d))
        direction_dir = os.path.join(dir_input, 'instrument_tables', d)
        file_list = glob.glob(os.path.join(direction_dir, '*'))
        inst_table_dir = os.path.join(dir_output, 'results', 'facetselfcal', d)
        for i, f in enumerate(file_list):
            log.info('  Copying instrument table file {0} of {1}...'.format(i+1, len(file_list)))
            copy(f, inst_table_dir, clobber)

        log.info('Unarchiving pipeline state files for direction {}...'.format(d))
        direction_dir = os.path.join(dir_input, 'pipeline_state', d)
        file_list = glob.glob(os.path.join(direction_dir, '*', '*'))
        mapfile_dir_list = []
        for i, f in enumerate(file_list):
            subdir = f.split('/')[-2]
            mapfile_dir = os.path.join(dir_output, 'results', subdir, d, 'mapfiles')
            mapfile_dir_list.append(mapfile_dir)
            copy(f, mapfile_dir, clobber)
        for mapfile_dir in set(mapfile_dir_list):
            update_state(mapfile_dir)

    log.info('Unarchiving complete.')
