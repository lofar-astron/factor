import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Changes the directory of files in the input mapfile to new_dir

    Parameters
    ----------
    mapfile_in : str
        Name of the input mapfile to be expanded. (E.g. with the skymodels for the
        different groups.)
    new_dir : str
        Name of the local directory to sync to
    make_tempdir : bool
        Make temporary directory inside new_dir?
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Output datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    new_dir = kwargs['new_dir']
    make_tempdir = True
    if 'make_tempdir' in kwargs:
        make_tempdir =  string2bool(kwargs['make_tempdir'])
    if make_tempdir:
        new_dir = tempfile.mkdtemp(dir=new_dir)
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_in.iterator = DataMap.SkipIterator
    map_out = DataMap([])
    for item in map_in:
        file_out = os.path.join(new_dir, os.path.dirname(item.file))
        map_out.data.append(DataProduct(item.host, file_out, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result


def string2bool(instring):
    if not isinstance(instring, basestring):
        raise ValueError('string2bool: Input is not a basic string!')
    if instring.upper() == 'TRUE' or instring == '1':
        return True
    elif instring.upper() == 'FALSE' or instring == '0':
        return False
    else:
        raise ValueError('string2bool: Cannot convert string "'+instring+'" to boolean!')
