import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Copies a mapfile

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap to copy
    hosts : str
        List of hosts/nodes. May be given as a list or as a string as a string
        (e.g., '[host1, host2]'
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    if 'hosts' in kwargs:
        if type(kwargs['hosts']) is str:
            hosts = kwargs['hosts'].strip('[]').split(',')
            hosts = [h.strip() for h in hosts]
        for i in range(len(datamap)-len(hosts)):
            hosts.append(hosts[i])
    else:
        hosts = None
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map = DataMap([])
    datamap = DataMap.load(mapfile_in)

    for i, item in enumerate(datamap):
        if hosts is None:
            host = item.host
        else:
            host = hosts[i]
        map.data.append(DataProduct(host, item.file, item.skip))

    fileid = os.path.join(mapfile_dir, filename)
    map.save(fileid)
    result = {'mapfile': fileid}

    return result
