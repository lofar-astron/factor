"""
Module that holds all cluster-related functions
"""
import os
import logging
import sys


log = logging.getLogger('cluster')


def make_pbs_clusterdesc(node_local_disk='/tmp'):
    """
    Make a cluster description file from the PBS_NODEFILE
    """
    nodes = []
    try:
        filename = os.environ['PBS_NODEFILE']
    except KeyError:
        self.log.error('PBS_NODEFILE not found. You must have a reservation to '
            'use clusterdesc = PBS.')
        sys.exit(1)

    with open(filename,'r') as file:
        for line in file:
            node_name = line.split()[0]
            if node_name not in nodes:
                nodes.append(node_name)

    lines = ['# Clusterdesc file to do parallel processing with PBS / torque\n\n']
    lines.append('ClusterName = PBS\n\n')
    lines.append('# Compute nodes\n')
    lines.append('Compute.Nodes = [{0}]\n'.format(', '.join(nodes)))
    lines.append('Compute.LocalDisks = [{0}]\n'.format(node_local_disk))

    clusterdesc_file = 'factor.clusterdesc'
    with open(clusterdesc_file, 'wb') as file:
        file.writelines(lines)

    return clusterdesc_file


def get_compute_nodes(clusterdesc_file):
    """
    Read a cluster description file and return list of nodes
    """
    from lofarpipe.support import clusterdesc

    cluster = clusterdesc.ClusterDesc(clusterdesc_file)
    return sorted(clusterdesc.get_compute_nodes(cluster))
