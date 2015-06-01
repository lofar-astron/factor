import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile for skymodel files given as a list

    Parameters
    ----------
    skymodel_files : str
        List of skymodel files (one per MS file) as a string (e.g.,
        '[s1.skymodel, s2.skymodel]'
    mapfile_in : str
        Filename of datamap containing MS files for input skymodels
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Skymodel datamap filename

    """
    skymodels = kwargs['skymodel_files'].strip('[]').split(',')
    skymodels = [s.strip() for s in skymodels]
    mapfile_in = kwargs['mapfile_in']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_out = DataMap([])

    map_in.iterator = DataMap.SkipIterator
    for i, item in enumerate(map_in):
        map_out.data.append(DataProduct(item.host, skymodels[i], False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
