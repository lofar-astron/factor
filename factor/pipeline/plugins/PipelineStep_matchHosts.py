import os
import glob
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Matchs the hosts in one datamap with those in another

    Parameters
    ----------
    mapfile_in : str, optional
        Filename of datamap to adjust
    mapfile_to_match : str, optional
        Filename of datamap to match

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_to_match = kwargs['mapfile_to_match']

    map_in = DataMap.load(mapfile_in)
    map_in.iterator = DataMap.SkipIterator
    map_to_match = DataMap.load(mapfile_to_match)
    map_to_match.iterator = DataMap.SkipIterator

    hosts_to_match = []
    for item in map_to_match:
        hosts_to_match.append(item.host)

    for item, host in zip(map_in, hosts_to_match):
        item.host = host

    map_in.save(mapfile_in)
