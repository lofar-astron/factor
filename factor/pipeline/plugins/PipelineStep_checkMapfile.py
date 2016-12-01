import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Checks a "check" mapfile for values of 'None' and, if found, changes the
    input mapfile "file" to "empty".

    Note: the check and input mapfiles must have the same length

    Parameters
    ----------
    mapfile_in : str
        Name of the input mapfile from which to select files.
    mapfile_check : str
        Name of the mapfile to check for None
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Output datamap filename

    """
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    inmap = DataMap.load(kwargs['mapfile_in'])
    checkmap = DataMap.load(kwargs['mapfile_check'])

    if len(inmap) != len(checkmap):
        raise ValueError('Input and check mapfiles must have the same length')

    map_out = DataMap([])
    for checkitem, item in zip(checkmap, inmap):
        if checkitem.file.lower() == 'none':
            map_out.data.append(DataProduct(item.host, 'empty', item.skip))
        else:
            map_out.append(item)

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result

