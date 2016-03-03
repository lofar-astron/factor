import os
import glob
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Updates the hosts in an input datamap

    Parameters
    ----------
    mapfile_in : str, optional
        Filename of datamap
    mapfile_dir: str, optional
        Directory containing mapfiles. All mapfiles in this directory will be
        updated
    hosts : str
        List of hosts/nodes. May be given as a list or as a string
        (e.g., '[host1, host2]'

    """
    if 'mapfile_dir' in kwargs:
        mapfiles_in = glob.glob(os.path.join(kwargs['mapfile_dir'], '*'))
    else:
        mapfiles_in = [kwargs['mapfile_in']]

    if len(mapfiles_in) == 0:
        return

    if type(kwargs['hosts']) is str:
        hosts = kwargs['hosts'].strip('[]').split(',')
        hosts = [h.strip() for h in hosts]

    for mapfile_in in mapfiles_in:
        try:
            map = DataMap.load(mapfile_in)
            for i in range(len(map)-len(hosts)):
                hosts.append(hosts[i])

            for item, host in zip(map, hosts):
                item.host = host

            map.save(mapfile_in)
        except:
            print('File {} does not appear to be a mapfile. Skipping it.'.format(mapfile_in))
