"""
Some functions called by multiple operations
"""
import logging
import os


def copy_column(ms_to, inputcol, outputcol, ms_from=None):
    """
    Copies one column to another

    Parameters
    ----------
    ms_to : str
        MS file receiving copy
    inputcol : str
        Column name to copy from
    outputcol : str
        Column name to copy to
    ms_from : str, optional
        MS file to copy from. If None, the column is copied internally

    """
    import pyrap.tables as pt
    import numpy as np

    t_to = pt.table(ms_to, readonly=False, ack=False)
    t0_to = np.min(t_to.getcol('TIME'))
    t1_to = np.max(t_to.getcol('TIME'))

    if ms_from is not None:
        t_from = pt.table(ms_from, ack=False)
        if outputcol not in t_to.colnames():
            cdesc = t_from.getcoldesc(inputcol)
            cdesc['name'] = outputcol
            t_to.addcols(cdesc)
        t0_from = np.min(t_from.getcol('TIME'))
        t1_from = np.max(t_from.getcol('TIME'))

        if t0_from >= t0_to and t1_from <= t1_to:
            # From-table is of equal or shorter length
            tinx_to = t_to.index(['TIME', "ANTENNA1", "ANTENNA2"])
            rownrs_to = tinx_to.rownrs({'TIME':t0_from, 'ANTENNA1':0, 'ANTENNA2':0},
                {'TIME':t1_from, 'ANTENNA1':500, 'ANTENNA2':500}, lowerincl=True,
                upperincl=True)
            tinx_from = t_from.index(['TIME', "ANTENNA1", "ANTENNA2"])
            rownrs_from = tinx_from.rownrs({'TIME':t0_from, 'ANTENNA1':0, 'ANTENNA2':0},
                {'TIME':t1_from, 'ANTENNA1':500, 'ANTENNA2':500}, lowerincl=True,
                upperincl=True)
        else:
            # From-table is longer
            tinx_to = t_to.index(['TIME', "ANTENNA1", "ANTENNA2"])
            rownrs_to = tinx_to.rownrs({'TIME':t0_to, 'ANTENNA1':0, 'ANTENNA2':0},
                {'TIME':t1_to, 'ANTENNA1':500, 'ANTENNA2':500}, lowerincl=True,
                upperincl=True)
            tinx_from = t_from.index(['TIME', "ANTENNA1", "ANTENNA2"])
            rownrs_from = tinx_from.rownrs({'TIME':t0_to, 'ANTENNA1':0, 'ANTENNA2':0},
                {'TIME':t1_to, 'ANTENNA1':500, 'ANTENNA2':500}, lowerincl=True,
                upperincl=True)

        data_to = t_to.getcol(outputcol)
        data_from = t_from.getcol(inputcol)
        data_to[rownrs_to] = data_from[rownrs_from]
    else:
        data_to = t_to.getcol(inputcol)
        if outputcol not in t_to.colnames():
            cdesc = t_to.getcoldesc(inputcol)
            cdesc['name'] = outputcol
            t_to.addcols(cdesc)

    # Write the changes
    t_to.putcol(outputcol, data_to)
    t_to.close()


def copy_column_freq(mslist, ms_from, inputcol, outputcol):
    """
    Copies column to bands

    Parameters
    ----------
    mslist : list
        MS files receiving copy
    ms_from : str
        MS file to copy from.
    inputcol : str
        Column name to copy from
    outputcol : str
        Column name to copy to

    """
    import sys
    import pyrap.tables as pt
    import numpy

    datain = pt.table(ms_from)
    data = datain.getcol(inputcol, nrow=1)
    numberofchans = numpy.int(numpy.shape(data)[1])
    chanperms = numberofchans/numpy.int(len(mslist))

    for ms_id, ms in enumerate(mslist):
        if os.path.isdir(ms):
            data = datain.getcolslice(inputcol, [chanperms*ms_id,0], [(chanperms*(ms_id+1))-1,3])
            dataout = pt.table(ms, readonly=False)
            dataout.putcol(outputcol, data)
            dataout.flush()
            dataout.close()


def make_chunks(dataset, blockl, op_parset, prefix=None, direction=None,
    columns=None, outdir=None, clobber=False):
    """
    Split dataset into time chunks of length chunksize time slots

    Parameters
    ----------
    dataset : str
        Name of MS file to split
    blockl : int
        Number of time slots per chunk
    op_parset : dict
        Parset of parent operation
    prefix : str
        A prefix for the name
    direction : Direction object or str, optional
        A direction name
    columns : str, optional
        List of column names to chunk. If None, all columns are chunked.
    outdir : str, optional
        Absolute path to output directory. If None, the directory of the
        parent is used
    clobber : bool, optional
        If True, existing files are overwritten

    """
    from factor.lib.chunk import Chunk
    import numpy as np
    import pyrap.tables as pt

    if blockl < 1:
        blockl = 1

    # Get time per sample and number of samples
    t = pt.table(dataset, readonly=True, ack=False)
    for t2 in t.iter(["ANTENNA1","ANTENNA2"]):
        if (t2.getcell('ANTENNA1',0)) < (t2.getcell('ANTENNA2',0)):
            timepersample = t2[1]['TIME']-t2[0]['TIME'] # sec
            nsamples = t2.nrows()
            break
    t.close()

    nchunks = int(np.ceil((np.float(nsamples) / np.float(blockl))))
    tlen = timepersample * np.float(blockl) / 3600. # length of block in hours
    tobs = timepersample * nsamples / 3600.0 # length of obs in hours

    # Set up the chunks
    chunk_list = []
    for c in range(nchunks):
        chunk_obj = Chunk(op_parset, dataset, c, prefix=prefix,
            direction=direction, outdir=outdir)
        chunk_obj.t0 = tlen * float(chunk_obj.index) # hours
        if c == 0:
            chunk_obj.t0 -= 0.1 # make sure first chunk gets first slot
        chunk_obj.t1 = np.float(chunk_obj.t0) + tlen # hours
        if c == nchunks-1 and chunk_obj.t1 < tobs:
            chunk_obj.t1 = tobs + 0.1 # make sure last chunk gets all that remains
        chunk_obj.start_delay = 0.0
        chunk_list.append(chunk_obj)
        split_ms(dataset, chunk_obj.file, chunk_obj.t0, chunk_obj.t1,
            columns=columns, clobber=clobber)

    return chunk_list


def split_ms(msin, msout, start_out, end_out, columns=None, clobber=True):
    """
    Splits an MS between start and end times in hours relative to first time

    Parameters
    ----------
    msin : str
        Name of MS file to split
    msout : str
        Name of output MS file
    start_out : float
        Start time in hours relative to first time
    end_out : float
        End time in hours relative to first time
    columns : str, optional
        List of column names to split. If None, all columns are split.
    clobber : bool, optional
        If True, existing files are overwritten

    """
    import pyrap.tables as pt
    import os

    if os.path.exists(msout):
        if clobber:
            os.system('rm -rf {0}'.format(msout))
        else:
            return

    t = pt.table(msin, ack=False)
    starttime = t[0]['TIME']

    if columns is not None:
        colnames = ', '.join(columns)
    else:
        colnames = ''

    t1 = t.query('TIME >= ' + str(starttime+start_out*3600) + ' && '
      'TIME < ' + str(starttime+end_out*3600), sortlist='TIME,ANTENNA1,ANTENNA2',
      columns=colnames)

    t1.copy(msout, True)
    t1.close()
    t.close()


def merge_chunks(chunk_files, prefix=None, virtual=True, clobber=False):
    """
    Merges chunks

    Parameters
    ----------
    chunk_files : list
        List of MS files to merge
    prefix : str, optional
        String to prepend to output file
    virtual : bool. optional
        If True, perform a virtual concat
    clobber : bool, optional
        If True, existing files are overwritten

    """
    import pyrap.tables as pt
    import os

    msout = None
    if prefix is not None:
        pstr = prefix + '_'
    else:
        pstr = ''
    for m in chunk_files:
        if '-chunk_0' in m:
            rstr = pstr + 'allchunks'
            msout = rstr.join(m.split('chunk_0'))
            break
    if msout is None:
        msout = chunk_files[0]

    if os.path.exists(msout):
        if clobber:
            os.system('rm -rf {0}'.format(msout))
        else:
            return msout

    if virtual:
        pt.msutil.msconcat(chunk_files, msout)
    else:
        t = pt.table(chunk_files, ack=False)
        t.sort('TIME').copy(msout, deep=True)
        t.close()

    return msout


def concatenate_chunks(chunk_files, prefix=None clobber=False):
    """
    Does a virtual oncatenation of chunks

    Parameters
    ----------
    chunk_files : list
        List of MS files to merge
    prefix: str, optional
        String to prepend to output file
    clobber : bool, optional
        If True, existing files are overwritten

    """
    msout = None
    if prefix is not None:
        pstr = prefix + '_'
    else:
        pstr = ''
    for m in chunk_files:
        if '-chunk_0' in m:
            rstr = pstr + 'allchunks'
            msout = rstr.join(m.split('chunk_0'))
            break
    if msout is None:
        msout = chunk_files[0]

    return concat_msname


def merge_chunk_parmdbs(inparmdbs, prefix='merged', clobber=False):
    """Merges parmdbs"""
    import os
    import lofar.parmdb
    import pyrap.tables as pt

    root_dir = inparmdbs[0].split('chunks')[0]
    outparmdb = '{0}/{1}_instrument'.format(root_dir, prefix)
    if os.path.exists(outparmdb):
        if clobber:
            os.system('rm -rf {0}'.format(outparmdb))
        else:
            return outparmdb

    os.system('cp -r {0} {1}'.format(inparmdbs[0], outparmdb))
    pdb_concat = lofar.parmdb.parmdb(outparmdb)

    for parmname in pdb_concat.getNames():
        pdb_concat.deleteValues(parmname)

    for inparmdb in inparmdbs:
        pdb = lofar.parmdb.parmdb(inparmdb)
        for parmname in pdb.getNames():
            v = pdb.getValuesGrid(parmname)
            pdb_concat.addValues(v)
    pdb_concat.flush()

    return outparmdb


def merge_parmdbs(parmdb_p, pre_apply_parmdb, parmdb_a, ms, clobber=True):
    """Merges amp+phase parmdbs"""
    import lofar.parmdb
    import pyrap.tables as pt
    import numpy as np

    parmdb_out = parmdb_p.split('phases2')[0] + 'final' + parmdb_p.split('phases2')[1]
    if os.path.exists(parmdb_out):
        if clobber:
            os.system('rm -rf {0}'.format(parmdb_out))
        else:
            return parmdb_out
    pdb_out = lofar.parmdb.parmdb(parmdb_out, create=True)

    pol_list = ['0:0', '1:1']
    gain = 'Gain'
    anttab = pt.table(ms + '::ANTENNA', ack=False)
    antenna_list = anttab.getcol('NAME')
    anttab.close()

    # Copy over the CommonScalar phases and TEC
    pdb_p = lofar.parmdb.parmdb(parmdb_p)
    for parmname in pdb_p.getNames():
        parms = pdb_p.getValuesGrid(parmname)
        pdb_out.addValues(parms)
        pdb_out.flush()

    # Get amplitude solutions
    pdb_pre = lofar.parmdb.parmdb(pre_apply_parmdb)
    pdb_a = lofar.parmdb.parmdb(parmdb_a)
    parms_pre = pdb_pre.getValuesGrid("*")
    parms_a = pdb_a.getValuesGrid("*")

    # Get array sizes and initialize values using first antenna (all
    # antennas should be the same)
    parmname = 'Gain:0:0:Real:' + antenna_list[0]
    N_times, N_freqs = parms_a[parmname]['values'].shape
    times = parms_a[parmname]['times'].copy()
    timewidths = parms_a[parmname]['timewidths'].copy()
    freqs = parms_a[parmname]['freqs'].copy()
    freqwidths = parms_a[parmname]['freqwidths'].copy()
    parms = {}
    v = {}
    v['times'] = times
    v['timewidths'] = timewidths
    v['freqs'] = freqs
    v['freqwidths'] = freqwidths

    # Multiply gains
    for pol in pol_list:
        for antenna in antenna_list:
            real1 = np.copy(parms_pre[gain+':'+pol+':Real:'+antenna]['values'])
            real2 = np.copy(parms_a[gain+':' +pol+':Real:'+antenna]['values'])
            imag1 = np.copy(parms_pre[gain+':'+pol+':Imag:'+antenna]['values'])
            imag2 = np.copy(parms_a[gain+':'+pol+':Imag:'+antenna]['values'])

            G1 = real1 + 1j * imag1
            G2 = real2 + 1j * imag2
            Gnew = G1 * G2
            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)

            parmname = gain + ':' + pol + ':Imag:' + antenna
            v['values'] = np.copy(np.imag(Gnew))
            parms[parmname] = v.copy()

            parmname = gain + ':' + pol + ':Real:' + antenna
            v['values'] = np.copy(np.real(Gnew))
            parms[parmname] = v.copy()

    pdb_out.addValues(parms)
    pdb_out.flush()

    return parmdb_out

