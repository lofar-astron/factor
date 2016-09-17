#! /usr/bin/env python
"""
Script to sort a list of MSs by into frequency groups by time-stamp
"""
import casacore.tables as pt
import sys, os
import numpy as np
import uuid
from lofarpipe.support.data_map import DataMap, DataProduct


def main(ms_input, filename=None, mapfile_dir=None, numSB=-1, enforce_numSB=True,
    hosts=None, NDPPPfill=True, target_path=None, stepname=None, nband_pad=0,
    make_dummy_files=False, skip_flagged_groups=True):
    """
    Check a list of MS files for missing frequencies

    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile
    filename: str
        Name of output mapfile
    mapfile_dir : str
        Directory for output mapfile
    numSB : int, optional
        How many files should go into one frequency group. Values <= 0 mean put
        all files of the same time-step into one group.
        default = -1
    enforce_numSB : bool, optional
        If True and numSB > 0, then add flagged dummy data to ensure that the
        last block has exactly numSB files. If False, then the last block can
        have fewer files (as long as there are no gaps in frequency)
    hosts : list or str
        List of hostnames or string with list of hostnames
    NDPPPfill : bool, optional
        Add dummy file-names for missing frequencies, so that NDPPP can
        fill the data with flagged dummy data.
        default = True
    target_path : str, optional
        Change the path of the "groups" files to this. (I.e. write output files
        into this directory with the subsequent NDPPP call.)
        default = keep path of input files
    stepname : str, optional
        Add this step-name into the file-names of the output files
    nband_pad : int, optional
        Add this number of bands of dummy data to the high-frequency end
        of the list
    make_dummy_files : bool, optional
        If True, make MS files for all dummy data
    skip_flagged_groups : bool, optional
        If True, groups that are missing have their skip flag set to True. If
        False, these groups are filled with dummy data and their skip flag set
        to False

    Returns
    -------
    result : dict
        Dict with the name of the generated mapfile

    """
    if not filename or not mapfile_dir:
        raise ValueError('sort_times_into_freqGroups: filename and mapfile_dir are needed!')

    # convert input to needed types
    ms_list = input2strlist(ms_input)
    NDPPPfill = input2bool(NDPPPfill)
    numSB = int(numSB)
    nband_pad = int(nband_pad)
    enforce_numSB = input2bool(enforce_numSB)
    make_dummy_files = input2bool(make_dummy_files)

    if type(hosts) is str:
        hosts = [h.strip(' \'\"') for h in hosts.strip('[]').split(',')]
    if not hosts:
        hosts = ['localhost']
    numhosts = len(hosts)
    print "sort_times_into_freqGroups: Working on",len(ms_list),"files"

    dirname = os.path.dirname(ms_list[0])

    time_groups = {}
    # sort by time
    for i, ms in enumerate(ms_list):
        # use the slower but more reliable way:
        obstable = pt.table(ms, ack=False)
        timestamp = int(round(np.min(obstable.getcol('TIME'))))
        #obstable = pt.table(ms+'::OBSERVATION', ack=False)
        #timestamp = int(round(obstable.col('TIME_RANGE')[0][0]))
        obstable.close()
        if timestamp in time_groups:
            time_groups[timestamp]['files'].append(ms)
        else:
            time_groups[timestamp] = {'files': [ ms ], 'basename' : os.path.splitext(ms)[0] }
    print "sort_times_into_freqGroups: found",len(time_groups),"time-groups"

    # sort time-groups by frequency
    timestamps = time_groups.keys()
    timestamps.sort()   # not needed now, but later
    first = True
    for time in timestamps:
        freqs = []
        for ms in time_groups[time]['files']:
            # Get the frequency info
            sw = pt.table(ms+'::SPECTRAL_WINDOW', ack=False)
            freq = sw.col('REF_FREQUENCY')[0]
            if first:
                freq_width = sw.col('TOTAL_BANDWIDTH')[0]
                maxfreq = freq
                minfreq = freq
                first = False
            else:
                assert freq_width == sw.col('TOTAL_BANDWIDTH')[0]
                maxfreq = max(maxfreq,freq)
                minfreq = min(minfreq,freq)
            freqs.append(freq)
            sw.close()
        time_groups[time]['freq_names'] = zip(freqs,time_groups[time]['files'])
        time_groups[time]['freq_names'].sort(key=lambda pair: pair[0])
        #time_groups[time]['files'] = [name for (freq,name) in freq_names]
        #time_groups[time]['freqs'] = [freq for (freq,name) in freq_names]
    print "sort_times_into_freqGroups: Collected the frequencies for the time-groups"

    #the new output map
    filemap = MultiDataMap()
    groupmap = DataMap()
    maxfreq = maxfreq+freq_width/2.
    minfreq = minfreq-freq_width/2.
    numFiles = round((maxfreq-minfreq)/freq_width)
    if numSB > 0:
        ngroups = int(np.ceil(numFiles/numSB))
    else:
        ngroups = 1
        numSB = int(numFiles)
    hostID = 0
    for time in timestamps:
        (freq,fname) = time_groups[time]['freq_names'].pop(0)
        nbands = 0
        all_group_files = []
        for fgroup in range(ngroups):
            files = []
            skip_this = True
            for fIdx in range(numSB):
                if freq > (fIdx+fgroup*numSB+1)*freq_width+minfreq:
                    files.append('dummy.ms')
                else:
                    files.append(fname)
                    skip_this = False
                    if len(time_groups[time]['freq_names'])>0:
                        (freq,fname) = time_groups[time]['freq_names'].pop(0)
                    else:
                        # Set freq to high value to pad the rest of the group
                        # with dummy data
                        (freq,fname) = (1e12,'This_shouldn\'t_show_up')
                        if not enforce_numSB:
                            # Don't pad the rest of this group with dummy data
                            break

            if fgroup == ngroups-1:
                # Append dummy data to last frequency group only
                for i in range(nband_pad):
                    files.append('dummy.ms')
            if not skip_this:
                nbands += len(files)

            if make_dummy_files:
                for i, ms in enumerate(files):
                    if ms == 'dummy.ms':
                        # Replace dummy.ms in files list with new filename
                        files[i] = os.path.join(dirname, '{0}_{1}.ms'.format(
                            os.path.splitext(ms)[0], uuid.uuid4().urn.split('-')[-1]))

            if not skip_flagged_groups:
                # Don't set skip flag to True, even if group is missing all files
                if not make_dummy_files:
                    raise ValueError('skip_flagged_groups cannot be False if make_dummy_files is also False')
                else:
                    skip_this = False

            filemap.append(MultiDataProduct(hosts[hostID%numhosts], files, skip_this))
            groupname = time_groups[time]['basename']+'_%Xt_%dg.ms'%(time,fgroup)
            if type(stepname) is str:
                groupname += stepname
            if type(target_path) is str:
                groupname = os.path.join(target_path,os.path.basename(groupname))
            groupmap.append(DataProduct(hosts[hostID%numhosts],groupname, skip_this))
            hostID += 1
            all_group_files.extend(files)

        assert freq==1e12

        if make_dummy_files:
            # Find at least one existing ms for this timestamp
            ms_exists = None
            for ms in all_group_files:
                if os.path.exists(ms):
                    ms_exists = ms
                    sw = pt.table('{}::SPECTRAL_WINDOW'.format(ms))
                    ms_exists_ref_freq = sw.getcol('REF_FREQUENCY')[0]
                    sw.close()
                    break

            for i, ms in enumerate(all_group_files):
                if 'dummy' in ms:
                    pt.tableutil.tablecopy(ms_exists, ms)

                    # Alter SPECTRAL_WINDOW subtable as appropriate to fill gap
                    sw = pt.table('{}::SPECTRAL_WINDOW'.format(ms), readonly=False)
                    tot_bandwidth = sw.getcol('TOTAL_BANDWIDTH')[0]
                    if i > 0:
                        sw_low = pt.table('{}::SPECTRAL_WINDOW'.format(all_group_files[i-1]))
                        ref_freq = sw_low.getcol('REF_FREQUENCY') + tot_bandwidth
                        sw_low.close()
                    else:
                        for j in range(1, len(files)-1):
                            if os.path.exists(files[j]):
                                sw_high = pt.table('{}::SPECTRAL_WINDOW'.format(all_group_files[j]))
                                ref_freq = sw_high.getcol('REF_FREQUENCY') - tot_bandwidth * j
                                sw_high.close()
                                break
                    chan_freq = sw.getcol('CHAN_FREQ') - ms_exists_ref_freq + ref_freq
                    sw.putcol('REF_FREQUENCY', ref_freq)
                    sw.putcol('CHAN_FREQ', chan_freq)
                    sw.close()

                    # Flag all data
                    t = pt.table(ms, readonly=False)
                    t.putcol('FLAG_ROW', np.ones(len(t), dtype=bool))
                    f = t.getcol('FLAG')
                    t.putcol('FLAG', np.ones(f.shape, dtype=bool))
                    t.close()

    filemapname = os.path.join(mapfile_dir, filename)
    filemap.save(filemapname)
    groupmapname = os.path.join(mapfile_dir, filename+'_groups')
    groupmap.save(groupmapname)
    result = {'mapfile': filemapname, 'groupmapfile': groupmapname, 'nbands': nbands}
    return result

def input2bool(invar):
    if isinstance(invar, bool):
        return invar
    elif isinstance(invar, str):
        if invar.upper() == 'TRUE' or invar == '1':
            return True
        elif invar.upper() == 'FALSE' or invar == '0':
            return False
        else:
            raise ValueError('input2bool: Cannot convert string "'+invar+'" to boolean!')
    elif isinstance(invar, int) or isinstance(invar, float):
        return bool(invar)
    else:
        raise TypeError('input2bool: Unsupported data type:'+str(type(invar)))

def input2strlist(invar):
    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
        else:
            map_in = DataMap.load(invar)
            map_in.iterator = DataMap.SkipIterator
            str_list = []
            for fname in map_in:
                if fname.startswith('[') and fname.endswith(']'):
                    for f in fname.strip('[]').split(','):
                        str_list.append(f.strip(' \'\"'))
                else:
                    str_list.append(fname.strip(' \'\"'))
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return str_list

class MultiDataProduct(DataProduct):
    """
    Class representing multiple files in a DataProduct.
    """
    def __init__(self, host=None, file=None, skip=True):
        super(MultiDataProduct, self).__init__(host, file, skip)
        if not file:
            self.file = list()
        else:
            self._set_file(file)

    def __repr__(self):
        """Represent an instance as a Python dict"""
        return (
            "{{'host': '{0}', 'file': '{1}', 'skip': {2}}}".format(self.host,
                '[{}]'.format(','.join(self.file)), str(self.skip))
        )

    def __str__(self):
        """Print an instance as 'host:[filelist]'"""
        return ':'.join((self.host, str(self.file)))

    def _set_file(self, data):
        try:
            # Try parsing as a list
            if isinstance(data, list):
                self.file = data
            if isinstance(data, DataProduct):
                self._from_dataproduct(data)
            if isinstance(data, DataMap):
                self._from_datamap(data)

        except TypeError:
            raise DataProduct("No known method to set a filelist from %s" % str(file))

    def _from_dataproduct(self, prod):
        print 'setting filelist from DataProduct'
        self.host = prod.host
        self.file = prod.file
        self.skip = prod.skip

    def _from_datamap(self, inmap):
        print 'setting filelist from DataMap'
        filelist = {}
        for item in inmap:
            if not item.host in filelist:
                filelist[item.host] = []
            filelist[item.host].append(item.file)
        self.file = filelist['i am']

    def append(self, item):
        self.file.append(item)


class MultiDataMap(DataMap):
    """
    Class representing a specialization of data-map, a collection of data
    products located on the same node, skippable as a set and individually
    """
    def __init__(self, data=list(), iterator=iter):
        super(MultiDataMap, self).__init__(data, iterator)

    @classmethod
    def load(cls, filename):
        """Load a data map from file `filename`. Return a DataMap instance."""
        with open(filename) as f:
            datamap = eval(f.read())
            file_entry = datamap['file']
            if file_entry.startswith('[') and file_entry.endswith(']'):
                file_list = [e.strip(' \'\"') for e in file_entry.strip('[]').split(',')]
                datamap = [{'host': datamap['host'], 'file': e, 'skip': datamap['skip']}
                    for e in file_list]
            return cls(datamap)

    @DataMap.data.setter
    def data(self, data):
        if isinstance(data, DataMap):
            mdpdict = {}
            data.iterator = DataMap.SkipIterator
            for item in data:
                if not item.host in mdpdict:
                    mdpdict[item.host] = []
                mdpdict[item.host].append(item.file)
            mdplist = []
            for k, v in mdpdict.iteritems():
                mdplist.append(MultiDataProduct(k, v, False))
            self._set_data(mdplist, dtype=MultiDataProduct)
        elif isinstance(data, MultiDataProduct):
            self._set_data(data, dtype=MultiDataProduct)
        elif not data:
            pass
        else:
            self._set_data(data, dtype=MultiDataProduct)

    def split_list(self, number):
        mdplist = []
        for item in self.data:
            for i in xrange(0, len(item.file), number):
                chunk = item.file[i:i+number]
                mdplist.append(MultiDataProduct(item.host, chunk, item.skip))
        self._set_data(mdplist)



if __name__ == '__main__':
    import optparse
    import glob
    import random

    opt = optparse.OptionParser(usage='%prog [options] <MSPattern> \n')
    opt.add_option('-v', '--verbose', help='Go Vebose! (default=False)', action='store_true', default=False)
    opt.add_option('-r', '--randomize', help='Randomize order of the input files. (default=False)', action='store_true', default=False)
    opt.add_option('-d', '--decimate', help='Remove every 10th file (after randomization if that is done). (default=False)', action='store_true', default=False)
    opt.add_option('-n', '--numbands', help='Number of how many files should be grouped together in frequency. (default=all files in one group)', type='int', default=-1)
    opt.add_option('-f', '--filename', help='Name for the mapfiles to write. (default=\"test.mapfile\")', type='string', default='test.mapfile')

    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 1:
        opt.print_help()
        sys.exit()

    # first argument: pattern for measurement-sets
    inMSs = glob.glob(args[0])
    if options.randomize:
        random.shuffle(inMSs)
    if options.decimate:
        for i in range((len(inMSs)-1),-1,-10):
            inMSs.pop(i)

    ergdict = main(inMSs, options.filename, '.', numSB=options.numbands, hosts=None, NDPPPfill=True)

    groupmap = DataMap.load(ergdict['groupmapfile'])
    filemap = MultiDataMap.load(ergdict['mapfile'])
    print "len(groupmap) : %d , len(filemap) : %d " % (len(groupmap),len(filemap))
    if len(groupmap) != len(filemap):
        print "groupmap and filemap have different length!"
        sys.exit(1)
    for i in xrange(len(groupmap)):
        print "Group \"%s\" has %d entries."%(groupmap[i].file,len(filemap[i].file))
