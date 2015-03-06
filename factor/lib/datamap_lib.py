"""
Datamap functions
"""
import os
from lofarpipe.support.data_map import DataMap, DataProduct


def write_mapfile(data_list, op_name, action_name=None, prefix=None,
    direction=None, index=None, working_dir='.', flag_list=None,
    use_abs_path=True, host_list=None):
    """
    Returns a datamap for the input data list

    Parameters
    ----------
    data_list : list of str
        List of files
    op_name : str
        Parent operation name
    action_name : str
        Parent action name
    prefix : str
        A prefix for the name
    direction : Direction object or str, optional
        A direction name
    index : int, optional
        An index for the datamap
    working_dir : str, optional
        Factor working directory
    flag_list : list of bools, optional
        List of datamap flags
    use_abs_path : bool, optional
        If True, make sure all files use absolute paths
    host_list : list, optional
        List of datamap host names

    """
    if use_abs_path:
        data_list = [os.path.abspath(f) for f in data_list]
    basename = make_mapfile_basename(prefix=prefix, direction=direction,
        index=index)
    mapfile_dir = os.path.join(working_dir, 'datamaps', op_name)
    if action_name is not None:
        mapfile_dir = os.path.join(mapfile_dir, action_name)
    if direction is not None:
        mapfile_dir = os.path.join(mapfile_dir, direction.name)
    if not os.path.exists(mapfile_dir):
        os.makedirs(mapfile_dir)
    mapfile = os.path.join(mapfile_dir, '{0}.datamap'.format(basename))

    datamap = DataMap([])
    if flag_list is None:
        flag_list = [False] * len(data_list)
    if host_list is None:
        host_list = ['localhost'] * len(data_list)
    for i in range(len(data_list)-len(host_list)):
        host_list.append(host_list[i])

    for data, flag, host in zip(data_list, flag_list, host_list):
        datamap.data.append(DataProduct(host, data, flag))
        datamap.save(mapfile)

    return mapfile


def read_mapfile(mapfile):
    """
    Returns a list of files and hosts in the input map file
    """
    indata = DataMap.load(mapfile)
    files = []
    hosts = []
    for item in indata:
        files.append(item.file)
        hosts.append(item.host)

    return files, hosts


def read_mapfile_flags(mapfile):
    """
    Returns a list of skip flags in the input map file

    Parameters
    ----------
    mapfile : str
        File name of Datamap file to read

    """
    indata = DataMap.load(mapfile)
    flags = []
    for item in indata:
        flags.append(item.skip)

    return flags


def set_mapfile_flags(mapfile, flag_list):
    """
    Sets flages in a datamap

    Parameters
    ----------
    mapfile : str
        File name of Datamap file to alter
    flag_list : list of bools
        List of flags to set

    """
    datamap = DataMap.load(mapfile)
    if len(datamap) != len(flag_list):
        print('Length of flag list does not match length of mapfile')
        return

    for item, flag in zip(datamap.data, flag_list):
        item.skip = flag

    datamap.save(mapfile)


def make_mapfile_basename(prefix=None, direction=None, index=None):
    """
    Returns a standard name pattern for datamap files

    Parameters
    ----------
    prefix : str
        A prefix for the name
    direction : Direction object, optional
        A direction name
    index : int, optional
        An index for the datamap

    """
    if direction is not None:
        dirtxt = '_{0}'.format(direction.name)
    else:
        dirtxt = ''
    if index is not None:
        indtxt = '-{0}'.format(index)
    else:
        indtxt = ''

    return '{0}{1}{2}'.format(prefix, dirtxt, indtxt)
