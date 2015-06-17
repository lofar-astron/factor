import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Removes files in a mapfile

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap to copy

    Returns
    -------
    result : dict
        New datamap filename

    """
    mapfile_in = kwargs['mapfile_in']

    datamap = DataMap.load(mapfile_in)

    for i, item in enumerate(datamap):
        if os.path.exists(item.file):
            os.system('rm -rf {0}'.format(item.file))
