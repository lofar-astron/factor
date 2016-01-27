import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by splitting a single input mapfile item into many items

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing single item
    split_str : str
        String to use to split
    hosts : list or str
        List of hosts/nodes. May be given as a list or as a string (e.g.,
        '[host1, host2]'
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    if 'split_str' not in kwargs:
        split_str = '||'
    if type(kwargs['hosts']) is str:
        hosts = kwargs['hosts'].strip('[]').split(',')
        hosts = [h.strip() for h in hosts]
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    files = mapfile_in[0].file
    files = files.split(split_str)
    files = [f.strip() for f in files]

    map_out = DataMap([])

    for i in range(len(files)-len(hosts)):
        hosts.append(hosts[i])

    map_out = DataMap([])
    for h, f in zip(hosts, files):
        map_out.data.append(DataProduct(h, f, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
