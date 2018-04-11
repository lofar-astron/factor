"""
Module that holds all compute-cluster-related functions
"""
import os
import logging
import sys
import numpy as np
from collections import Counter
import factor._logging
import re

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


def expand_part(s):
    """Expand a part (e.g. "x[1-2]y[1-3][1-3]") (no outer level commas).

    Note: Adapted from git://www.nsc.liu.se/~kent/python-hostlist.git
    """

    # Base case: the empty part expand to the singleton list of ""
    if s == "":
        return [""]

    # Split into:
    # 1) prefix string (may be empty)
    # 2) rangelist in brackets (may be missing)
    # 3) the rest

    m = re.match(r'([^,\[]*)(\[[^\]]*\])?(.*)', s)
    (prefix, rangelist, rest) = m.group(1,2,3)

    # Expand the rest first (here is where we recurse!)
    rest_expanded = expand_part(rest)

    # Expand our own part
    if not rangelist:
        # If there is no rangelist, our own contribution is the prefix only
        us_expanded = [prefix]
    else:
        # Otherwise expand the rangelist (adding the prefix before)
        us_expanded = expand_rangelist(prefix, rangelist[1:-1])

    return [us_part + rest_part
            for us_part in us_expanded
            for rest_part in rest_expanded]


def expand_range(prefix, range_):
    """ Expand a range (e.g. 1-10 or 14), putting a prefix before.

    Note: Adapted from git://www.nsc.liu.se/~kent/python-hostlist.git
    """

    # Check for a single number first
    m = re.match(r'^[0-9]+$', range_)
    if m:
        return ["%s%s" % (prefix, range_)]

    # Otherwise split low-high
    m = re.match(r'^([0-9]+)-([0-9]+)$', range_)

    (s_low, s_high) = m.group(1,2)
    low = int(s_low)
    high = int(s_high)
    width = len(s_low)

    results = []
    for i in xrange(low, high+1):
        results.append("%s%0*d" % (prefix, width, i))
    return results


def expand_rangelist(prefix, rangelist):
    """ Expand a rangelist (e.g. "1-10,14"), putting a prefix before.

    Note: Adapted from git://www.nsc.liu.se/~kent/python-hostlist.git
    """

    # Split at commas and expand each range separately
    results = []
    for range_ in rangelist.split(","):
        results.extend(expand_range(prefix, range_))
    return results


def expand_hostlist(hostlist, allow_duplicates=False, sort=False):
    """Expand a hostlist expression string to a Python list.

    Example: expand_hostlist("n[9-11],d[01-02]") ==>
             ['n9', 'n10', 'n11', 'd01', 'd02']

    Unless allow_duplicates is true, duplicates will be purged
    from the results. If sort is true, the output will be sorted.

    Note: Adapted from git://www.nsc.liu.se/~kent/python-hostlist.git
    """

    results = []
    bracket_level = 0
    part = ""

    for c in hostlist + ",":
        if c == "," and bracket_level == 0:
            # Comma at top level, split!
            if part: results.extend(expand_part(part))
            part = ""
            bad_part = False
        else:
            part += c

        if c == "[": bracket_level += 1
        elif c == "]": bracket_level -= 1

    seen = set()
    results_nodup = []
    for e in results:
        if e not in seen:
            results_nodup.append(e)
            seen.add(e)
    return results_nodup


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
        hostlist = os.environ['SLURM_JOB_NODELIST']
    except KeyError:
        log.error('SLURM_JOB_NODELIST not found. You must have a reservation to '
            'use clusterdesc = SLURM.')
        sys.exit(1)

    nodes = expand_hostlist(hostlist)

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
                   'image2fits_executable': ['image2fits'],
                   'h5collector_executable': ['H5parm_collector.py']}
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
