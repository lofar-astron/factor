import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by filtering input mapfile items into one item (the middle
    one)

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
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
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_out = DataMap([])

    map_in.iterator = DataMap.SkipIterator
    files = [item.file for item in map_in]
    hosts = [item.host for item in map_in]
    if 'index' in kwargs:
        index = int(kwargs['index'])
    else:
        index = len(files)/2
    map_out.data.append(DataProduct(hosts[index], files[index], False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
