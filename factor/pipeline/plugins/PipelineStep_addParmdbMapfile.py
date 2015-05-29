import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile for parmdb files in an input MS datamap

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
    parmdb_name : str
        Name of parmdb files inside the MS files
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
    parmdb_name = kwargs['parmdb_name']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_out = DataMap([])

    map_in.iterator = DataMap.SkipIterator
    for item in map_in:
        map_out.data.append(DataProduct(item.host, os.path.join(item.file,
            parmdb_name), False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
