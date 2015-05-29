import os
from lofarpipe.support.data_map import DataMap


def plugin_main(args, **kwargs):
    """
    Copies a mapfile

    Parameters
    ----------
    mapfile : str
        Filename of MS datamap to copy
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New datamap filename

    """
    mapfile = kwargs['mapfile']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    datamap = DataMap.load(mapfile)

    fileid = os.path.join(mapfile_dir, filename)
    datamap.save(fileid)
    result = {'mapfile': fileid}

    return result
