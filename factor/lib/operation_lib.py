"""
Some functions called by multiple operations
"""
import logging
import os


def copy_column(ms, inputcol, outputcol, ms_from=None,):
    """
    Copies one column to another

    Parameters
    ----------
    ms : str
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

    t = pt.table(ms, readonly=False, ack=False)

    t0 = t[0]['TIME']
    t1 = t[-1]['TIME']

    if ms_from is not None:
        tf = pt.table(ms_from, readonly=False, ack=False)
        cd = tf.getcoldesc(inputcol)
        cd['name'] = outputcol
        try:
            t.addcols(cd)
        except RuntimeError:
            # Column already exists
            pass
        t0f = tf[0]['TIME']
        t1f = tf[-1]['TIME']
        if t0f >= t0 and t1f <= t1:
            # From-table is of equal or shorter length, so get indices for
            # to-table
            startrow = np.where(t.getcol('TIME') >= t0f)[0][0]
            nrow = np.where(t.getcol('TIME') <= t1f)[0][-1] - startrow + 1
            data = tf.getcol(inputcol)
            t.putcol(inputcol, data, startrow=startrow, nrow=nrow)
        else:
            # From-table is longer, so get indices for from-table
            startrow = np.where(tf.getcol('TIME') >= t0)[0][0]
            nrow = np.where(tf.getcol('TIME') < t1)[0][-1] - startrow
            data = tf.getcol(inputcol, startrow=startrow, nrow=nrow)
            t.putcol(outputcol, data)
    else:
        data = t.getcol(inputcol)
        cd = t.getcoldesc(inputcol)
        cd['name'] = outputcol
        try:
            t.addcols(cd)
        except:
            pass
        t.putcol(outputcol, data)

    # Write the changes
    t.flush()
    t.close()


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


def merge_chunks(chunk_files, prefix=None, clobber=False):
    """
    Merges chunks

    Parameters
    ----------
    chunk_files : list
        List of MS files to merge
    prefix: str, optional
        String to prepend to output file
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

    t = pt.table(chunk_files, ack=False)
    t.sort('TIME').copy(msout, deep = True)
    t.close()
    return msout
