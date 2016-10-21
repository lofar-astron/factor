import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Appends a string to filenames in a mapfile

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap to append to
    append : str
        String to append
    append_index : bool
        If True, append a unique index to each file
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

    if 'append_index' in kwargs:
        append_index = kwargs['append_index']
        if type(append_index) is str:
            if append_index.lower() == 'true':
                append_index = True
            else:
                append_index = False
    else:
        append_index = False

    append_str = kwargs['append']
    if append_str == 'None':
        append_str = ''
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_out = DataMap([])
    map_in = DataMap.load(mapfile_in)

    for i, item in enumerate(map_in):
        if append_index:
            map_out.data.append(DataProduct(item.host, item.file+append_str+'_{}'.format(i), item.skip))
        else:
            map_out.data.append(DataProduct(item.host, item.file+append_str, item.skip))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
