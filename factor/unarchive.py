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

log = logging.getLogger('factor:unarchive')


def update_state(mapfile_dir):
    """
    Updates the paths in state files

    Parameters
    ----------
    dir_input : str
        Path to files to update

    """
    pipe_dir = os.path.dirname(mapfile_dir)
    file_list = glob.glob(os.path.join(mapfile_dir, '*'))
    for f in file_list:
        # Try to load file as a pipeline mapfile
        try:
            map = DataMap.load(f)
            for item in map:
                if '/' in item.file:
                    item.file = os.path.join(pipe_dir, os.path.basename(item.file))
            map.save(f)
        except:
            # Try to load file as a pickled file and look for paths inside
            try:
                fp = open(f, 'r')
                d = pickle.load(fp)
                fp.close()
                for k, v in d.iteritems():
                    if type(v) is str:
                        if 'mapfiles' in v:
                            parts = v.split('mapfiles')
                            d[k] = os.path.join(pipe_dir, 'mapfiles', parts[1])
                        if 'skymodel_dirindep' in v:
                            filename = os.path.basename(v)
                            bandname = d['name']
                            d[k] = os.path.join(pipe_dir.split('state')[0],
                                'sky_models', bandname, filename)
                    elif type(v) is list:
                        for i, l in enumerate(v):
                            if '/chunks/' in l:
                                parts = l.split('/chunks/')
                                v[i] = os.path.join(pipe_dir, 'chunks', parts[1])
                        d[k] = v

                fp = open(f, 'w')
                pickle.dump(d, fp)
                fp.close()
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
        for i, f in enumerate(file_list):
            subdir = f.split('/')[-2]
            mapfile_dir = os.path.join(dir_output, 'results', subdir, d, 'mapfiles')
            copy(f, mapfile_dir, clobber)
        update_state(mapfile_dir)

    log.info('Unarchiving complete.')
