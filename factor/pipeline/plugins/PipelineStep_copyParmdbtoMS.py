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

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    result = {}
    ms_datamap = kwargs[ms_mapfile]
    parmdb_datamap = kwargs[parmdb_mapfile]
    for ms, parmdb in zip(ms_datamap, parmdb_datamap):
        from_file = parmdb.file
        to_file = os.path.join(ms.file, 'instrument')
        if os.path.exists(to_file):
            os.system('rm -rf {0}'.format(to_file))
        os.system('cp -r {0} {1}'.format(from_file, to_file))
        parmdb.file = to_file

    out_mapfile = os.path.splitext(parmdb_mapfile)[0] + '_copyparmdbtoms.' + os.path.splitext(parmdb_mapfile)[1]
    parmdb_datamap.save(out_mapfile)
    result['mapfile'] = out_mapfile
    return result
