import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Copies a files in a mapfile

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap to copy
    prefix : str
        Prefix to add to original file names to get new ones
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
    prefix = kwargs['prefix']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map = DataMap([])
    datamap = DataMap.load(mapfile_in)

    for i, item in enumerate(datamap):
        new_file = os.path.join(os.path.dirname(item.file), prefix+os.path.basename(item.file))
        if os.path.exists(new_file):
            os.system('rm -rf {0}'.format(new_file))
        os.system('cp -r {0} {1}'.format(item.file, new_file))

        map.data.append(DataProduct(item.host, new_file, item.skip))

    fileid = os.path.join(mapfile_dir, filename)
    map.save(fileid)
    result = {'mapfile': fileid}

    return result
