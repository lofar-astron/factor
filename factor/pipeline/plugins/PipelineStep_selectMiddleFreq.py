import os
import pyrap.tables as pt
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile with only the MSs at the middle Frequency
    Quite a bit of a hack for a plugin, but right now I don't care.

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_in.iterator = DataMap.SkipIterator
    map_out = DataMap()

    #sort into frequency groups
    freq_groups = {}
    hosts = []
    for item in map_in:
        # Get the frequency info
        sw = pt.table(item.file+'::SPECTRAL_WINDOW', ack=False)
        freq = int(sw.col('REF_FREQUENCY')[0])
        sw.close()
        if freq in freq_groups:
            freq_groups[freq].append(item.file)
        else:
            freq_groups[freq] = [item.file]
        if not item.host in hosts:
            hosts.append(item.host)
    # find maximum number of files per frequency-group
    maxfiles = max([len(group) for group in freq_groups.values()])
    # find the center-frequency
    freqs =  freq_groups.keys()
    freqs.sort()
    midfreq = freqs[len(freqs)/2]
    # make sure that chosen frequncy has maxfiles entries
    while len(freq_groups[midfreq]) < maxfiles:
        freqs.remove(midfreq)
        midfreq = freqs[len(freqs)/2]
    # extend the hosts-list
    for i in range(maxfiles-len(hosts)):
        hosts.append(hosts[i])
    # fill the output-map
    for (host,fname) in zip(hosts,freq_groups[midfreq]):
        map_out.append(DataProduct(host, fname, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
