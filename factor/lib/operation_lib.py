"""
Some functions called by multiple operations
"""
import logging
import os


def copy_column(ms, inputcol, outputcol, ms_from=None):
    """
    Copies one column to another
    """
    import sys
    import pyrap.tables as pt
    import os
    import numpy

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

