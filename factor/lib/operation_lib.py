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
    outdir=None, clobber=False):
    """
    Split dataset into time chunks of length blockl time slots

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
    tlen = timepersample * np.float(blockl) / 3600.0 # length of block in hours
    tobs = timepersample * nsamples / 3600.0 # length of obs in hours

    # Set up the chunks
    chunk_list = []
    for c in range(nchunks):
        chunk_obj = Chunk(op_parset, dataset, c, prefix=prefix,
            direction=direction, outdir=outdir)
        chunk_obj.t0 = tlen * np.float(chunk_obj.index) # hours
        chunk_obj.t1 = chunk_obj.t0 + tlen # hours
        if c == 0:
            chunk_obj.t0 = -0.1 # make sure first chunk gets first slot
        if c == nchunks-1 and chunk_obj.t1 < tobs:
            chunk_obj.t1 = tobs + 0.1 # make sure last chunk gets all that remains
        chunk_obj.start_delay = 0.0
        chunk_list.append(chunk_obj)
        split_ms(dataset, chunk_obj.file, chunk_obj.t0, chunk_obj.t1,
            clobber=clobber)

    return chunk_list


def split_ms(msin, msout, start_out, end_out, clobber=True):
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

    t1 = t.query('TIME >= ' + str(starttime+start_out*3600.0) + ' && '
      'TIME < ' + str(starttime+end_out*3600.0), sortlist='TIME,ANTENNA1,ANTENNA2')

    t1.copy(msout, True)
    t1.close()
    t.close()


def merge_chunks(chunk_files, prefix=None, virtual=True, clobber=False):
    """
    Merges chunks

    Parameters
    ----------
    chunk_files : list
        List of files to merge
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
        msout = pstr + chunk_files[0] + '_merged'

    if os.path.exists(msout):
        if clobber:
            os.system('rm -rf {0}'.format(msout))
        else:
            return msout

    if virtual:
        # TODO: put in order by time?
        pt.msutil.msconcat(chunk_files, msout)
    else:
        t = pt.table(chunk_files, ack=False)
        t.sort('TIME').copy(msout, deep=True)
        t.close()

    return msout


def merge_chunk_parmdbs(inparmdbs, prefix='merged', clobber=True):
    """
    Merges chunk parmdbs into a single parmdb

    Parameters
    ----------
    inparmdbs : list
        List of input parmdb file names
    prefix : str, optional
        String to prefix to output parmdb file name
    clobber : bool, optional
        If True, overwrite existing output file

    """
    import lofar.parmdb

    root_dir = inparmdbs[0].split('chunks')[0]
    outparmdb = '{0}/{1}_instrument'.format(root_dir, prefix)
    if os.path.exists(outparmdb):
        if clobber:
            os.system('rm -rf {0}'.format(outparmdb))
        else:
            return outparmdb

    os.system('cp -r {0} {1}'.format(inparmdbs[0], outparmdb))

    if len(inparmdbs) > 1:
        pdb_concat = lofar.parmdb.parmdb(outparmdb)
        for inparmdb in inparmdbs[1:]:
            pdb = lofar.parmdb.parmdb(inparmdb)
            for parmname in pdb.getNames():
                v = pdb.getValuesGrid(parmname)
                pdb_concat.addValues(v.copy())
        pdb_concat.flush()

    return outparmdb


def check_selfcal(image_prev, image_final, max_rms, max_ratio):
    """
    Checks that selfcal is improving
    """
    return False


def merge_parmdbs(parmdb_p, parmdb_a, prefix=None, clobber=True):
    """
    Merges facet selfcal parmdbs into a parmdb for a single band

    Parameters
    ----------
    parmdb_p : str
        File name of CommonScalarPhase and TEC parmdb
    parmdb_a : str
        File name of Gain parmdb. The nearset match in frequency to that of the
        input band will be used
    prefix : str, optional
        Prefix for name of output file
    clobber : bool, optional
        If True, overwrite existing output file

    """
    import lofar.parmdb

    # Initialize output parmdb
    if prefix is None:
        prefix = 'merged'
    parmdb_out = os.path.join(os.path.dirname(parmdb_p),
        '{0}_amp_phase_final_instrument'.format(prefix))
    if os.path.exists(parmdb_out):
        if clobber:
            os.system('rm -rf {0}'.format(parmdb_out))
        else:
            return parmdb_out
    pdb_out = lofar.parmdb.parmdb(parmdb_out, create=True)

    # Copy over the CommonScalar phases and TEC
    pdb_p = lofar.parmdb.parmdb(parmdb_p)
    for parmname in pdb_p.getNames():
        parms = pdb_p.getValuesGrid(parmname)
        pdb_out.addValues(parms)

    # Copy over the Gains
    pdb_a = lofar.parmdb.parmdb(parmdb_a)
    for parmname in pdb_a.getNames():
        parms = pdb_a.getValuesGrid(parmname)
        pdb_out.addValues(parms)

    # Write values
    pdb_out.flush()

    return parmdb_out


def verify_subtract(image_pre, image_post, res_val, imager):
    """
    Check quantities in residual images
    """
    import numpy
    import pyrap.images as pim

    if imager.lower() == 'wsclean':
        image_pre = image_pre + '-image.fits'
        image_post = image_post + '-image.fits'
    else:
        # Not implemented yet
        return True

    imgpre = pim.image(image_pre)
    pixelspre = numpy.copy(imgpre.getdata())
    maxvalpre = numpy.copy(numpy.max(pixelspre))

    img = pim.image(image_post)
    pixels = numpy.copy(img.getdata())
    maxval = numpy.copy(numpy.max(pixels))

    if (maxval > res_val) or ((maxval*0.95) > maxvalpre) :
        logging.info('WARNING RESIDUAL TOO LARGE')
        logging.info('Max = {0}'.format(maxval))
        logging.info('Previous max = {0}'.format(maxvalpre))
        return False
    else:
        return True
