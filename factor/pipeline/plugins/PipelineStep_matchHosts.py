import os
import glob
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Matchs the hosts in one datamap with another

    Parameters
    ----------
    mapfile_in : str, optional
        Filename of datamap to adjust
    mapfile_to_match : str, optional
        Filename of datamap to match
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_to_match = kwargs['mapfile_to_match']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_in.iterator = DataMap.SkipIterator
    map_to_match = DataMap.load(mapfile_to_match)
    map_to_match.iterator = DataMap.SkipIterator
    map_out = DataMap([])

    hosts_to_match = []
    for item in map_to_match:
        hosts_to_match.append(item.host)

    for item, host in zip(map_in, hosts_to_match):
        map_out.data.append(DataProduct(host, item.file, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
