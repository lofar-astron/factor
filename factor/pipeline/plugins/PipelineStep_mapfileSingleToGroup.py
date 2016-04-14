import os
from lofarpipe.support.data_map import DataMap, MultiDataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Copies each entry of mapfile_in as often as the the length of the corresponding
    group into a new mapfile

    Parameters
    ----------
    mapfile_in : str
        Name of the input mapfile to be expanded. (E.g. with the skymodels for the
        different groups.)
    mapfile_groups : str
        Name of the multi-mapfile with the given groups. Number of groups need
        to be the same as the number of files in mapfile_in.
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
    groupmap = MultiDataMap.load(kwargs['mapfile_groups'])

    if len(inmap) != len(groupmap):
        raise ValueError('PipelineStep_mapfileSingleToGroup: length of {0} and {1} differ'.format(kwargs['mapfile_in'],kwargs['mapfile_groups']))

    map_out = DataMap([])
    inindex = 0
    for groupID in xrange(len(groupmap)):
        for fileID in xrange(len(groupmap[groupID].file)):
            map_out.data.append(DataProduct(inmap[groupID].host, inmap[groupID].file, (inmap[groupID].skip or groupmap[groupID].skip) ))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
