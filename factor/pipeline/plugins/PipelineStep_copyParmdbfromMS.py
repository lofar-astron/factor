import os
from lofarpipe.support.data_map import DataMap

def plugin_main(args, **kwargs):
    """
    Copies parmdb from MS/{{suffix}}instrument to new location

    Parameters
    ----------
    ms_mapfile : str
        Filename of MS datamap
    parmdb_name : str
        Name of parmdb files inside the MS files
    suffix : str
        Suffix to append to intput parmdb name to get output name
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    ms_mapfile = kwargs['ms_mapfile']
    parmdb_name = kwargs['parmdb_name']
    suffix = kwargs['suffix']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    ms_datamap = DataMap.load(ms_mapfile)

    for ms in ms_datamap:
        from_file = os.path.join(ms.file, parmdb_name)
        to_file = from_file + suffix
        if os.path.exists(to_file):
            os.system('rm -rf {0}'.format(to_file))
        os.system('cp -r {0} {1}'.format(from_file, to_file))
        ms.file = to_file

    fileid = os.path.join(mapfile_dir, filename)
    ms_datamap.save(fileid)
    result = {'mapfile': fileid}

    return result
