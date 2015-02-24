"""
Datamap functions
"""
import os
from lofarpipe.support.data_map import DataMap, DataProduct


def write_mapfile(data_list, op_name, action_name=None, prefix=None,
    direction=None, index=None, working_dir='.'):
    """
    Returns a datamap for the input data list
    """
    data_list = [os.path.abspath(f) for f in data_list]
    basename = make_mapfile_basename(action_name, prefix, direction, index)
    mapfile = os.path.join(working_dir, 'datamaps', op_name,
        '{0}.datamap'.format(basename))

    datamap = DataMap([])
    for data in data_list:
        datamap.data.append(DataProduct('localhost', data, False))
        datamap.save(mapfile)

    return mapfile


def read_mapfile(mapfile):
    """
    Returns a list of files in the input map file
    """
    indata = DataMap.load(mapfile)
    files = []
    for item in indata:
        files.append(item.file)

    return files


def make_mapfile_basename(action_name=None, prefix=None, direction=None, index=None):
    """
    Returns a standard name pattern for datamap files
    """
    if direction is not None:
        try:
            dirtxt = '_{0}'.format(direction.name)
        except:
            dirtxt = '_{0}'.format(direction)
    else:
        dirtxt = ''
    if index is not None:
        indtxt = '-{0}'.format(index)
    else:
        indtxt = ''

    if action_name is not None:
        return '{0}/{1}{2}{3}'.format(action_name, prefix, dirtxt, indtxt)
    else:
        return '{0}{1}{2}'.format(prefix, dirtxt, indtxt)
