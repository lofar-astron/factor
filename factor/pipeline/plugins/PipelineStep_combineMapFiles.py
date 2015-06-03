import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile of by combining items from input mapfiles into one item

    Parameters
    ----------
    mapfiles_in : list or str
        List of filenames of datamaps. May be given as a list or as a string (e.g.,
        '[datamap1, datamap2]'
    list_of_str : bool
        If True, combined item will be a list of strings (e.g., "['file1',
        'file2']"). Mainly used for casapy calls
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    if type(kwargs['mapfiles_in']) is str:
        mapfiles_in = kwargs['mapfiles_in'].strip('[]').split(',')
        mapfiles_in = [m.strip() for m in mapfiles_in]
    num_mapfiles = len(mapfiles_in)
    list_of_str = bool(kwargs['list_of_str'])
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_out = DataMap([])

    files_list = []
    for mapfile_in in mapfiles_in:
        map_in = DataMap.load(mapfile_in)
        map_in.iterator = DataMap.SkipIterator
        files_list.append([item.file for item in map_in])
        hosts = [item.host for item in map_in]

    new_list = []
    for i in range(len(files_list[0]):
        for j in range(num_mapfiles):
            new_list.append(files_list[j][i])

        if list_of_str:
            newlist = ['"{0}"'.format(f) for f in newlist]
        map_out.data.append(DataProduct(hosts[i], '[{0}]'.format(','.join(newlist)), False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
