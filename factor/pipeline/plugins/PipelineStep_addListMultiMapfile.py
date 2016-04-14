import os
from lofarpipe.support.data_map import MultiDataMap, MultiDataProduct


def plugin_main(args, **kwargs):
    """
    Makes a multi-mapfile for list of files

    Parameters
    ----------
    files : str
        Nested list of files given as a string. The sub-lists _must_ be separated by the string "], ["
        '[[file1a, file1b, file1c], [file2a, file2b]]'
    hosts : str
        List of hosts/nodes. May be given as a list or as a string (e.g.,
        '[host1, host2]'
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Output datamap filename

    """
    if type(kwargs['files']) is str:
        tmplist = kwargs['files'].strip('[]').split('], [')
        files = []
        for filestring in tmplist:
            fileslist = filestring.split(',')
            files.append([f.strip(' \'') for f in fileslist])
    else:
        print 'PipelineStep_addListMultiMapfile.py kwargs[\'files\'] is not a string!'
        raise ValueError('PipelineStep_addListMultiMapfile.py kwargs[\'files\'] is not a string!')
    if type(kwargs['hosts']) is str:
        hosts = kwargs['hosts'].strip('[]').split(',')
        hosts = [h.strip() for h in hosts]
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    for i in range(len(files)-len(hosts)):
        hosts.append(hosts[i])

    map_out = MultiDataMap([])
    for h, f in zip(hosts, files):
        map_out.data.append(MultiDataProduct(h, f, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result

