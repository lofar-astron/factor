import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Re-groups a simple (flat) mapfile into a multi-mapfile that has the same shape
    as a given multi-mapfile.

    Parameters
    ----------
    mapfile_in : str
        Name of the input mapfile to be re-grouped.
    mapfile_groups : str
        Name of the multi-mapfile with the given groups. Total number of files needs 
        to be the same as in mapfile_in. 
    check_basename : Bool (str) , optional
        Check if the basenames (see os.path.basename()) minus extension match
        default = True
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
    check_names = True
    if 'check_basename' in kwargs:
        check_names =  string2bool(kwargs['check_basename'])

    inmap = DataMap.load(kwargs['mapfile_in'])
    groupmap = MultiDataMap.load(kwargs['mapfile_groups'])


    map_out = MultiDataMap([])
    inindex = 0
    for group in groupmap:
        grouplist = []
        skip = False
        for filename in group.file:
            if check_names:
                refbase = os.path.splitext(os.path.basename(filename))[0]
                newbase = os.path.splitext(os.path.basename(inmap[inindex].file))[0]
                if refbase != newbase:
                    raise ValueError('PipelineStep_reGroupMapfile: basenames {0} and {1} differ'.format(refbase,newbase))
            grouplist.append(inmap[inindex].file)
            if inmap[inindex].skip:
                print 'PipelineStep_reGroupMapfile: Skipping full group for file:'+inmap[inindex].file
                skip = True
            inindex += 1            
        map_out.data.append(MultiDataProduct(group.host, grouplist, skip))
    assert inindex == len(inmap)

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
            "{{'host': '{0}', 'file': {1}, 'skip': {2}}}".format(self.host, self.file, str(self.skip))
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
            print 'HELP: ', data
            self._set_data(data, dtype=MultiDataProduct)

    def split_list(self, number):
        mdplist = []
        for item in self.data:
            for i in xrange(0, len(item.file), number):
                chunk = item.file[i:i+number]
                mdplist.append(MultiDataProduct(item.host, chunk, item.skip))
        self._set_data(mdplist)

