import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by compressing input mapfile items into one item

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile
    list_format : bool, optional
        If True, the compreseed item will use a Python list format (e.g.,
        '[file1, file2, ...]'. If False, it will be a space-separated list (e.g.,
        'file1 file2 ...'

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    if 'list_format' in kwargs:
        list_format = kwargs['list_format']
    else:
        list_format = True
    if type(list_format) is str:
        if list_format.lower() == 'true':
            list_format = True
        else:
            list_format = False

    map_in = DataMap.load(mapfile_in)
    map_out = DataMap([])

    map_in.iterator = DataMap.SkipIterator
    if list_format:
        newlist = '[{0}]'.format(','.join([item.file for item in map_in]))
    else:
        newlist = '{0}'.format(' '.join([item.file for item in map_in]))
    hosts = [item.host for item in map_in]
    map_out.data.append(DataProduct(hosts[0], newlist, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
