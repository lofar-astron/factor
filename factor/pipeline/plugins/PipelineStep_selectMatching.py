import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Selects those files from mapfile_in that have the same filename-base as the one in
    mapfile_reference.

    Parameters
    ----------
    mapfile_in : str
        Name of the input mapfile from which to select files.
    mapfile_reference : str
        Name of the reference mapfile
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
    refmap = DataMap.load(kwargs['mapfile_reference'])

    map_out = DataMap([])

    basenames = [ os.path.splitext(os.path.basename(item.file))[0] for item in inmap]
    for refitem in refmap:
        refbase = os.path.splitext(os.path.basename(refitem.file))[0]
        idx = basenames.index(refbase)
        map_out.append(inmap[idx])

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result

