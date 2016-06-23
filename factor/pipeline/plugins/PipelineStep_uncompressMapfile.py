import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by uncompressing input mapfile list item into separate items

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing list of MS files
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile
    hosts : str
        List of hosts/nodes. May be given as a list or as a string
        (e.g., '[host1, host2]'

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    if type(kwargs['hosts']) is str:
        hosts = kwargs['hosts'].strip('[]').split(',')
        hosts = [h.strip() for h in hosts]

    map_in = DataMap.load(mapfile_in)
    map_out = DataMap([])

    files = map_in[0].file.strip('[]').split(',')
    files = [f.strip() for f in files]
    for i in range(len(files)-len(hosts)):
        hosts.append(hosts[i])

    for file, host in zip(files, hosts):
        map_out.data.append(DataProduct(host, file, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
