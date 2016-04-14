import os
import numpy as np
from lofarpipe.support.data_map import DataMap, DataProduct, MultiDataMap, MultiDataProduct


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
    join_groups : int (str), optional
        If it is set, then join so many groups into one new group. (Gives fewer
        groups but more files per group than in mapfile_groups.)
        default = keep same grouping as in mapfile_groups
    join_max_files : int (str), optional
        If it is set, then try to join as many groups together before the number of
        files per group woud exceed "join_max_files". Similar to "join_groups", but
        the number of joind groups is not fixed but depends on the number of files
        per group. Mutaully exclusive with "join_groups"!
    rotate_groups : bool (str), optional
        If set to "True": form groups by taking n files from each group, with n given
        by join_groups or join_max_files. I.e. form groups of files from all bands
        instead of groups with the same bands. This will fully mangle the order of
        the time-steps but will guarantee at least one file per (frequency-)group
        in each output-group.
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Output datamap filename

    """
    if 'join_groups' in kwargs and 'join_max_files' in kwargs:
        raise ValueError("PipelineStep_reGroupMapfile: \"join_groups\" and \"join_max_files\" are mutually exclusive!")

    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    check_names = True
    if 'check_basename' in kwargs:
        check_names =  string2bool(kwargs['check_basename'])
    rotate_groups = False
    if 'rotate_groups' in kwargs:
        rotate_groups = string2bool(kwargs['rotate_groups'])

    inmap = DataMap.load(kwargs['mapfile_in'])
    groupmap = MultiDataMap.load(kwargs['mapfile_groups'])

    map_out = MultiDataMap([])
    inindex = 0
    minpergroup = len(inmap)
    for group in groupmap:
        grouplist = []
        skip = False
        for fname in group.file:
            if check_names:
                refbase = os.path.splitext(os.path.basename(fname))[0]
                newbase = os.path.splitext(os.path.basename(inmap[inindex].file))[0]
                if refbase != newbase:
                    raise ValueError('PipelineStep_reGroupMapfile: basenames {0} and {1} differ'.format(refbase,newbase))
            grouplist.append(inmap[inindex].file)
            if inmap[inindex].skip:
                print 'PipelineStep_reGroupMapfile: Skipping full group for file:'+inmap[inindex].file
                skip = True
            inindex += 1
        map_out.data.append(MultiDataProduct(group.host, grouplist, skip))
        minpergroup = min(minpergroup,len(grouplist))
    assert inindex == len(inmap)


    if not rotate_groups:
        if 'join_groups' in kwargs:
            groups_to_join =  int(kwargs['join_groups'])
            newmap = MultiDataMap([])
            for start_idx in xrange(0,len(map_out),groups_to_join):
                end_idx = min((start_idx+groups_to_join),len(map_out))
                grouplist = []
                for group in map_out[start_idx:end_idx]:
                    grouplist.extend(group.file)
                    if group.skip:
                        raise ValueError("PipelineStep_reGroupMapfile: Found group that should be skipped! "
                                         "(I.e. there is probably something wrong with your data!)")
                newmap.data.append(MultiDataProduct(map_out[start_idx].host, grouplist, False))
            map_out = newmap
        elif 'join_max_files' in kwargs:
            max_files =  int(kwargs['join_max_files'])
            newmap = MultiDataMap([])
            grouplist = map_out[0].file
            grouphost = map_out[0].host
            for gindex in xrange(1,len(map_out)):
                if map_out[gindex].skip:
                    raise ValueError("PipelineStep_reGroupMapfile: Found group that should be skipped! "
                                     "(I.e. there is probably something wrong with your data!)")
                if (len(grouplist)+len(map_out[gindex].file)) > max_files:
                    newmap.data.append(MultiDataProduct(grouphost, grouplist, False))
                    grouplist = map_out[gindex].file
                    grouphost = map_out[gindex].host
                else:
                    grouplist.extend(map_out[gindex].file)
            # add the final (partial?) group to the map
            newmap.data.append(MultiDataProduct(grouphost, grouplist, False))
            map_out = newmap
    else:
       if 'join_groups' in kwargs or 'join_max_files' in kwargs:
           if 'join_groups' in kwargs:
               num_groups_out = max(1, minpergroup / int(kwargs['join_groups']))
           else:
               groups_to_join = np.ceil(float(kwargs['join_max_files'])/len(map_out))
               num_groups_out = max(1,int(minpergroup/groups_to_join))
           outgroups = {}
           for i in xrange(num_groups_out):
               outgroups[i] = []
           hosts = []
           for group in map_out:
               hosts.append(group.host)
               groupidxs = np.array_split(np.arange(len(group.file)),num_groups_out)
               for i in xrange(num_groups_out):
                   for j in groupidxs[i]:
                       outgroups[i].append(group.file[j])
           newmap = MultiDataMap([])
           for i in xrange(num_groups_out):
               newmap.data.append(MultiDataProduct(hosts[i%len(hosts)], outgroups[i], False))
           map_out = newmap

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
