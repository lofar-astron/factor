"""
Module that holds all cluster-related functions
"""
import os
import logging
import sys
import factor._logging


log = logging.getLogger('factor.cluster')


def make_pbs_clusterdesc(node_local_disk='/tmp'):
    """
    Make a cluster description file from the PBS_NODEFILE
    """
    nodes = []
    try:
        filename = os.environ['PBS_NODEFILE']
    except KeyError:
        log.error('PBS_NODEFILE not found. You must have a reservation to '
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
    lines.append('Compute.Nodes = [{0}]\n'.format(', '.join(sorted(nodes))))
    lines.append('Compute.LocalDisks = [{0}]\n'.format(node_local_disk))

    clusterdesc_file = 'factor.clusterdesc'
    with open(clusterdesc_file, 'wb') as file:
        file.writelines(lines)
    log.debug('Using {0} node(s)'.format(len(nodes)))

    return clusterdesc_file


def get_compute_nodes(clusterdesc_file):
    """
    Read a cluster description file and return list of nodes
    """
    from lofarpipe.support import clusterdesc

    cluster = clusterdesc.ClusterDesc(clusterdesc_file)
    return sorted(clusterdesc.get_compute_nodes(cluster))


def find_executables(parset):
    """
    Finds paths to required executables
    """
    from distutils import spawn

    executables = {'casa_executable': ['casa', 'casapy'],
                   'wsclean_executable': ['wsclean'],
#                    'chgcentre_executable': ['chgcentre'],
                   'losoto_executable': ['losoto.py'],
                   'H5parm_importer_executable': ['H5parm_importer.py'],
                   'H5parm_exporter_executable': ['H5parm_exporter.py']}
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
