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

