"""
Some functions called by multiple operations
"""
import logging
import os


def copy_column(ms, inputcol, outputcol, ms_from=None):
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

    t = pt.table(ms, readonly=False, ack=False)
    if ms_from is not None:
        tf = pt.table(ms_from, readonly=False, ack=False)
        data = tf.getcol(inputcol)
        cd = tf.getcoldesc(inputcol)
    else:
        data = t.getcol(inputcol)
        cd = t.getcoldesc(inputcol)
    cd['name'] = outputcol
    try:
        t.addcols(cd)
    except:
        pass
    t.putcol(outputcol, data)
    t.flush()
    t.close()


def make_chunks(dataset, blockl, op_parset, prefix=None, direction=None, clobber=False):
    """
    Split ms into time chunks of length chunksize time slots
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
            direction=direction)
        chunk_obj.t0 = tlen * float(chunk_obj.index) # hours
        chunk_obj.t1 = np.float(chunk_obj.t0) + tlen # hours
        if c == nchunks-1 and chunk_obj.t1 < tobs:
            chunk_obj.t1 = tobs + 0.1 # make sure last chunk gets all that remains
        chunk_obj.start_delay = 0.0
        chunk_list.append(chunk_obj)
        split_ms(dataset, chunk_obj.file, chunk_obj.t0, chunk_obj.t1,
            clobber=clobber)

    return chunk_list


def split_ms(msin, msout, start_out, end_out, clobber=True):
    """Splits an MS between start and end times in hours relative to first time"""
    import pyrap.tables as pt
    import os

    if os.path.exists(msout):
        if clobber:
            os.system('rm -rf {0}'.format(msout))
        else:
            return

    t = pt.table(msin, ack=False)

    starttime = t[0]['TIME']
    t1 = t.query('TIME > ' + str(starttime+start_out*3600) + ' && '
      'TIME < ' + str(starttime+end_out*3600), sortlist='TIME,ANTENNA1,ANTENNA2')

    t1.copy(msout, True)
    t1.close()
    t.close()


def merge_chunks(chunk_files, prefix=None, clobber=False):
    """Merges chunks"""
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
