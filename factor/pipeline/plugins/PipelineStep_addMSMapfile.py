import os
import glob
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile for MS files in a folder

    Parameters
    ----------
    folder : str
        Directory containing MS files
    hosts : list or str
        List of hosts/nodes. May be given as a list or as a string (e.g.,
        '[host1, host2]'
    mapfile_dir : str
        Output directory for mapfile
    filename: str
        Name of mapfile

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    folder = kwargs['folder']
    if type(kwargs['hosts']) is str:
        hosts = kwargs['hosts'].strip('[]').split(',')
        hosts = [h.strip() for h in hosts]
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_out = DataMap([])
    measurements = glob.glob(os.path.join(folder, '*[MS|ms]'))
    measurements.sort()

    for i in range(len(measurements)-len(hosts)):
        hosts.append(hosts[i])

    for ms, host in zip(measurements, hosts):
        map_out.data.append(DataProduct(host, os.path.join(folder, ms), False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
