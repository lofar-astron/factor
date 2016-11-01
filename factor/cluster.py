"""
Module that holds all compute-cluster-related functions
"""
import os
import logging
import sys
import numpy as np
from collections import Counter
import factor._logging


log = logging.getLogger('factor:cluster')


def make_pbs_clusterdesc():
    """
    Make a cluster description file from the PBS_NODEFILE

    Returns
    -------
    clusterdesc_file
        Filename of resulting cluster description file

    """
    nodes = []
    try:
        filename = os.environ['PBS_NODEFILE']
    except KeyError:
        log.error('PBS_NODEFILE not found. You must have a reservation to '
            'use clusterdesc = PBS.')
        sys.exit(1)

    with open(filename, 'r') as file:
        for line in file:
            node_name = line.split()[0]
            if node_name not in nodes:
                nodes.append(node_name)

    lines = ['# Clusterdesc file to do parallel processing with PBS / torque\n\n']
    lines.append('ClusterName = PBS\n\n')
    lines.append('# Compute nodes\n')
    lines.append('Compute.Nodes = [{0}]\n'.format(', '.join(sorted(nodes))))

    clusterdesc_file = 'factor_pbs.clusterdesc'
    with open(clusterdesc_file, 'wb') as file:
        file.writelines(lines)
    log.info('Using {0} node(s)'.format(len(nodes)))

    return clusterdesc_file


def make_slurm_clusterdesc():
    """
    Make a cluster description file from the SLURM_JOB_NODELIST

    Returns
    -------
    clusterdesc_file
        Filename of resulting cluster description file

    """
    nodes = []
    try:
        filename = os.environ['SLURM_JOB_NODELIST']
    except KeyError:
        log.error('SLURM_JOB_NODELIST not found. You must have a reservation to '
            'use clusterdesc = SLURM.')
        sys.exit(1)

    with open(filename, 'r') as file:
        for line in file:
            node_name = line.split()[0]
            if node_name not in nodes:
                nodes.append(node_name)

    lines = ['# Clusterdesc file to do parallel processing with SLURM\n\n']
    lines.append('ClusterName = SLURM\n\n')
    lines.append('# Compute nodes\n')
    lines.append('Compute.Nodes = [{0}]\n'.format(', '.join(sorted(nodes))))

    clusterdesc_file = 'factor_slurm.clusterdesc'
    with open(clusterdesc_file, 'wb') as file:
        file.writelines(lines)
    log.info('Using {0} node(s)'.format(len(nodes)))

    return clusterdesc_file


def get_compute_nodes(clusterdesc_file):
    """
    Read a cluster description file and return list of nodes

    Parameters
    ----------
    clusterdesc_file : str
        Filename of cluster description file

    Returns
    -------
    result : list
        Sorted list of node names

    """
    from lofarpipe.support import clusterdesc

    cluster = clusterdesc.ClusterDesc(clusterdesc_file)
    return sorted(clusterdesc.get_compute_nodes(cluster))


def find_executables(parset):
    """
    Adds the paths to required executables to parset dict

    Parameters
    ----------
    parset : dict
        Parset dictionary

    """
    from distutils import spawn

    executables = {'genericpipeline_executable': ['genericpipeline.py'],
                   'wsclean_executable': ['wsclean'],
                   'image2fits_executable': ['image2fits']}
    for key, names in executables.iteritems():
        for name in names:
            path = spawn.find_executable(name)
            if path is not None:
                parset[key] = path
                break
        if path is None:
            log.error('The path to the {0} executable could not be determined. '
                'Please make sure it is in your PATH.'.format(name))
            sys.exit(1)
