import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Updates the hosts in an input datamap

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap
    hosts : str
        List of hosts/nodes. May be given as a list or as a string
        (e.g., '[host1, host2]'

    Returns
    -------
    result : dict
        Input datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    if type(kwargs['hosts']) is str:
        hosts = kwargs['hosts'].strip('[]').split(',')
        hosts = [h.strip() for h in hosts]

    map = DataMap.load(mapfile_in)
    for i in range(len(map)-len(hosts)):
        hosts.append(hosts[i])

    for item, host in zip(map, hosts):
        item.host = host

    map.save(mapfile_in)
    result = {'mapfile': mapfile_in}

    return result
