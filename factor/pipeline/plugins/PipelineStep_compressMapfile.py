import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by compressing input mapfile items into one item

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile
    list_format : bool, optional
        If True, the compreseed item will use a Python list format (e.g.,
        '[file1, file2, ...]'. If False, it will be a space-separated list (e.g.,
        'file1 file2 ...'

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    if 'list_format' in kwargs:
        list_format = kwargs['list_format']
    else:
        list_format = True
    if type(list_format) is str:
        if list_format.lower() == 'true':
            list_format = True
        else:
            list_format = False

    map_in = DataMap.load(mapfile_in)
    map_out = DataMap([])
    map_in.iterator = DataMap.SkipIterator

    # Try to detect missing bands by common naming conventions:
    # SB010 or B10
    # This assumes bands are ordered from low to high
    band_numbers = []
    file_list = []
    skip_check = False
    for i, item in enumerate(map_in):
        try:
            if 'SB' in item.file:
                band_indx = int(item.file.split('SB')[1][0:3])
            elif 'B' in item.file:
                band_indx = int(item.file.split('B')[1][0:2])
            if i == 0:
                start_indx = band_indx
            elif i == 1:
                indx_skip = band_indx - start_indx
            band_numbers.append(band_indx)
            file_list.append(item.file)
        except ValueError:
            file_list = [item.file for item in map_in]
            skip_check = True
            break
    if not skip_check:
        indx = 0
        for i in range(len(band_numbers)):
            if band_numbers[i] != (indx * indx_skip) + start_indx:
                file_list.insert(indx, 'dummy.ms')
                indx += 1
            indx += 1

    if list_format:
        newlist = '[{0}]'.format(','.join(file_list))
    else:
        newlist = '{0}'.format(' '.join(file_list))

    # Just assign host of first file to compressed file
    hosts = [item.host for item in map_in]
    map_out.data.append(DataProduct(hosts[0], newlist, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
