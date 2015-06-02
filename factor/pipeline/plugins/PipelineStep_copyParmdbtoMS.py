import os
from lofarpipe.support.data_map import DataMap

def plugin_main(args, **kwargs):
    """
    Copies parmdb to MS/instrument

    Parameters
    ----------
    ms_mapfile : str
        Filename of MS datamap
    parmdb_mapfile : str
        Filename of parmdb datamap

    """
    ms_mapfile = kwargs['ms_mapfile']
    parmdb_mapfile = kwargs['parmdb_mapfile']

    ms_datamap = DataMap.load(ms_mapfile)
    parmdb_datamap = DataMap.load(parmdb_mapfile)

    for ms, parmdb in zip(ms_datamap, parmdb_datamap):
        from_file = parmdb.file
        to_file = os.path.join(ms.file, 'instrument')
        if os.path.exists(to_file):
            os.system('rm -rf {0}'.format(to_file))
        os.system('cp -r {0} {1}'.format(from_file, to_file))
