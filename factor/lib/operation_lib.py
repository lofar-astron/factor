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


def merge_chunks(chunks, prefix=None, virtual=True, clobber=False):
    """
    Merges chunks

    Parameters
    ----------
    chunks : list
        List of Chunk objects to merge
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

    chunk_files = [chunk.file for chunk in chunks]
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
        # put in order by time
        pt.msutil.msconcat(chunk_files, msout)
    else:
        t = pt.table(chunk_files, ack=False)
        t.sort('TIME').copy(msout, deep=True)
        t.close()

    return msout


def merge_chunk_parmdbs(inparmdbs, prefix='merged', clobber=False):
    """
    Merges chunk parmdbs into a single parmdb
    """
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


def merge_parmdbs(parmdb_p, parmdb_a, parmdb_t, solint_p, solint_a, msfile,
    prefix=None, clobber=True):
    """
    Merges facet selfcal parmdbs into a parmdb for a single band

    Parameters
    ----------
    parmdb_p : str
        File name of CommonScalarPhase and TEC parmdb
    parmdb_a : str
        File name of Gain parmdb. The nearset match in frequency to that of the
        input band will be used
    parmdb_t : str
        File name of template parmdb
    solint_p : int
        Solution interval for parmdb_p
    solint_a : int
        Solution interval for parmdb_a
    msfile : str
        File name of band
    prefix : str
        Prefix to prepend to output file
    clobber : bool, optional
        If True, overwrite existing output file

    """
    import lofar.parmdb
    import pyrap.tables as pt
    import numpy as np

    # Initialize output parmdb
    parmdb_out = parmdb_p.split('phases2')[0] + 'final' + parmdb_p.split('phases2')[1]
    if os.path.exists(parmdb_out):
        if clobber:
            os.system('rm -rf {0}'.format(parmdb_out))
        else:
            return parmdb_out
    pdb_out = lofar.parmdb.parmdb(parmdb_out, create=True)

    # Open input parmdbs
    pdb_a = lofar.parmdb.parmdb(parmdb_a)
    pdb_p = lofar.parmdb.parmdb(parmdb_p)
    pdb_t = lofar.parmdb.parmdb(parmdb_t)

    # Get solutions
    parms_a = pdb_a.getValuesGrid("*")
    parms_p = pdb_p.getValuesGrid("*")
    parms_t = pdb_t.getValuesGrid("*")

    # Get antenna names
    pol_list = ['0:0', '1:1']
    anttab = pt.table(msfile + '::ANTENNA', ack=False)
    antenna_list = anttab.getcol('NAME')
    anttab.close()

    # Set time and frequency grid for output parmdb
    # The time grid can be set to that of parmdb_p (fast phase grid)
    # The freq grid must be set to that of parmdb_t (template grid)
    parmname = 'TEC:' + antenna_list[0]
    times = parms_p[parmname]['times'].copy()
    timewidths = parms_p[parmname]['timewidths'].copy()
    N_times_p, N_freqs_p = parms_p[parmname]['values'].shape

    if 'Gain' in pdb_t.getNames()[0]:
        parmname = 'Gain:0:0:Real:' + antenna_list[0]
    elif 'Phase' in pdb_t.getNames()[0]:
        parmname = 'Phase:0:0:' + antenna_list[0]
    freqs = parms_t[parmname]['freqs'].copy()
    freqwidths = parms_t[parmname]['freqwidths'].copy()
    N_times_t, N_freqs_t = parms_t[parmname]['values'].shape

    parmname = 'Gain:0:0:Real:' + antenna_list[0]
    N_times_a, N_freqs_a = parms_a[parmname]['values'].shape
    freqs_a = parms_a[parmname]['freqs'].copy()
    freq_ind = np.searchsorted(freqs_a, freqs)

    # Initialize parms and values dicts
    outparms_p = {}
    v_p = {}
    v_p['times'] = times
    v_p['timewidths'] = timewidths
    v_p['freqs'] = freqs
    v_p['freqwidths'] = freqwidths
    outparms_g = {}
    v_g = {}
    v_g['times'] = times
    v_g['timewidths'] = timewidths
    v_g['freqs'] = freqs
    v_g['freqwidths'] = freqwidths

    # Copy values
    for pol in pol_list:
        for antenna in antenna_list:
            # Copy gains
            v_g['values'] = np.zeros((N_times_p, N_freqs_t), dtype=np.double)

            parmname = 'Gain:' + pol + ':Imag:' + antenna
            imag = np.copy(parms_a[parmname]['values'])[:, freq_ind]
            imag_repeat = np.repeat(imag, solint_a/solint_p, axis=0)
            v_g['values'] = np.copy(imag_repeat[0:N_times_p])
            outparms_g[parmname] = v_g.copy()

            parmname = 'Gain:' + pol + ':Real:' + antenna
            real = np.copy(parms_a[parmname]['values'])[:, freq_ind]
            real_repeat = np.repeat(real, solint_a/solint_p, axis=0)
            v_g['values'] = np.copy(real_repeat[0:N_times_p])
            outparms_g[parmname] = v_g.copy()

            # Copy CommonScalar phases and TEC
            v_p['values'] = np.zeros((N_times_p, N_freqs_t), dtype=np.double)

            parmname = 'CommonScalarPhase:' + antenna
            phase = np.copy(parms_p[parmname]['values'])
            v_p['values'] = np.copy(phase)
            outparms_p[parmname] = v_p.copy()

            parmname = 'TEC:' + antenna
            phase = np.copy(parms_p[parmname]['values'])
            v_p['values'] = np.copy(phase)
            outparms_p[parmname] = v_p.copy()

    pdb_out.addValues(outparms_g)
    pdb_out.addValues(outparms_p)
    pdb_out.flush()

    return parmdb_out


def verify_subtract(image_pre, image_post, res_val):
    """
    Check quantities in residual images
    """
    import numpy
    import pyrap.images

    imgpre = pyrap.images.image(image_pre)
    pixelspre = numpy.copy(imgpre.getdata())
    maxvalpre = numpy.copy(numpy.max(pixelspre))

    img = pyrap.images.image(image_post)
    pixels = numpy.copy(img.getdata())
    maxval = numpy.copy(numpy.max(pixels))

    if (maxval > res_val) or ((maxval*0.95) > maxvalpre) :
        print 'WARNING RESIDUAL TOO LARGE, STOPPING', maxval, res_val
        print 'WARNING RESIDUAL TOO LARGE, STOPPING, previous max in image', maxvalpre
        return False
    else:
        return True
