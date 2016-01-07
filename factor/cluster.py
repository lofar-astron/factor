"""
Module that holds all compute-cluster-related functions
"""
import os
import logging
import sys
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

    clusterdesc_file = 'factor.clusterdesc'
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

    executables = {'genericpipeline_executable': ['genericpipeline2.py', 'genericpipeline.py'],
                   'casa_executable': ['casapy', 'casa'],
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


def divide_nodes(directions, node_list, ndir_per_node, nimg_per_node, ncpu_max,
    fmem_max, nbands):
    """
    Divide up nodes and cpus for parallel operations

    Parameters
    ----------
    directions: list
        List of Direction objects for which the nodes will be divided
    node_list : list
        List of node names to divide among the directions
    ndir_per_node : int
        Number of directions per node
    nimg_per_node : int
        Number of images per node
    ncpu_max : int
        Max number of CPUs per node
    fmem_max : float
        Max fraction of memory per node
    nbands : int
        Number of bands

    Returns
    -------
    directions : list
        List of Direction objects, with attributes modified

    """
    if len(directions) >= len(node_list):
        for i in range(len(directions)-len(node_list)):
            node_list.append(node_list[i])
        hosts = [[n] for n in node_list]
    else:
        parts = len(directions)
        hosts = [node_list[i*len(node_list)//parts:
            (i+1)*len(node_list)//parts] for i in range(parts)]

    # Find duplicates and divide up available cores
    h_flat = []
    for h in hosts:
        h_flat.extend(h)
    c = Counter(h_flat)
    for d, h in zip(directions, hosts):
        d.hosts = h
        if len(h) == 1:
            ndir_per_node = min(ndir_per_node, c[h[0]])
        else:
            ndir_per_node = 1
        d.nimg_per_node = nimg_per_node
        d.max_cpus_per_node =  max(1, int(round(ncpu_max / float(ndir_per_node))))
        d.max_cpus_per_img =  max(1, int(round(ncpu_max / float(nimg_per_node))))
        nchunks_per_node = max(1, int(round(float(d.nchunks) / len(d.hosts))))
        d.max_cpus_per_chunk = int(round(d.max_cpus_per_node / nchunks_per_node))
        d.max_cpus_per_band = max(1, int(round(d.max_cpus_per_node *
            len(d.hosts) / float(nbands))))
        d.max_percent_memory = fmem_max / float(ndir_per_node) / float(nimg_per_node) * 100.0
        d.save_state()

    return directions


def combine_nodes(directions, node_list, nimg_per_node, ncpu_max, fmem_max, nbands):
    """
    Conbine nodes and cpus for serial operations

    Parameters
    ----------
    directions: list
        List of Direction objects
    node_list : list
        List of node names
    nimg_per_node : int
        Number of images per node
    ncpu_max : int
        Max number of CPUs per node
    fmem_max : float
        Max fraction of memory per node
    nbands : int
        Number of bands

    Returns
    -------
    directions : list
        List of Direction objects, with attributes modified

    """
    # Give each direction all resources
    for d in directions:
        d.hosts = node_list
        d.nimg_per_node = nimg_per_node
        d.max_cpus_per_node = ncpu_max
        d.max_cpus_per_img =  max(1, int(round(ncpu_max / float(nimg_per_node))))
        nchunks_per_node = max(1, int(round(float(d.nchunks) / len(d.hosts))))
        d.max_cpus_per_chunk = max(1, int(round(d.max_cpus_per_node /
            float(nchunks_per_node))))
        d.max_cpus_per_band = max(1, int(round(ncpu_max * len(d.hosts) / float(nbands))))
        d.max_percent_memory = fmem_max / float(nimg_per_node) * 100.0
        d.save_state()

    return directions
