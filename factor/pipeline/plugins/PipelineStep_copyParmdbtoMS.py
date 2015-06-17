import os
from lofarpipe.support.data_map import DataMap

def plugin_main(args, **kwargs):
    """
    Copies a single parmdb to MS/instrument

    A single parmdb can be copied to multiple MS files

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
    if len(parmdb_datamap) == 1:
        parmdb_files = [parmdb_datamap[0].file] * len(ms_datamap)
    else:
        parmdb_files = [item.file for item in parmdb_datamap]

    for ms, parmdb_file in zip(ms_datamap, parmdb_files):
        to_file = os.path.join(ms.file, 'instrument')
        if os.path.exists(to_file):
            os.system('rm -rf {0}'.format(to_file))
        os.system('cp -r {0} {1}'.format(parmdb_file, to_file))
