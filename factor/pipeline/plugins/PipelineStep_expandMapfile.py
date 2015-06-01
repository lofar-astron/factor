import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by expanding input mapfile item

    Note: only input mapfiles with a single item are supported

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap
    number : str
        Number of items in new mapfile
    suffix : str
        Suffix to add to new item files. The index of each item is appended to
        the suffix (e.g., suffix of '_chunk' -> file_chunk0.ext, file_chunk1.ext,
        etc.)
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
    number = kwargs['number']
    suffix = kwargs['suffix']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_out = DataMap([])

    for item in map_in:
        host = item.host
        newlist = ['{0}{1}{2}.ms'.format(os.path.splitext(item.file)[0], suffix, c)
            for c in range(number)]
    for newfile in newlist:
        map_out.data.append(DataProduct(host, newfile, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
