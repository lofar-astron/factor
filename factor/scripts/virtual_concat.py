#!/usr/bin/env python
"""
Script to perform a virtual concatenation
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import os


def main(ms_files, outfile, clobber=True):
    """
    Performs a virtual concatenation

    Parameters
    ----------
    ms_files : list
        List of files to merge
    outfile : str
        Output file
    clobber : bool, optional
        If True, existing files are overwritten

    """
    if type(ms_files) is str:
        ms_files = [f.strip() for f in ms_files.strip('[]').split(',')]
    if type(clobber) is str:
        if clobber.lower() == 'true':
            clobber = True
        else:
            clobber = False
    if os.path.exists(outfile):
        if clobber:
            os.system('rm -rf {0}'.format(outfile))
        else:
            return

    # TODO: order by time?
    try:
        pt.msutil.msconcat(ms_files, outfile)
    except ImportError:
        if os.path.exists(outfile):
            os.system('rm -rf {0}'.format(outfile))
        msconcat(ms_files, outfile)


def msconcat(names, newname, concatTime=False):
    """Virtually concatenate multiple MeasurementSets

    NOTE: this function was taken from the LOFAR trunk and altered to bypass
    an import problem on CEP3

    Multiple MeasurementSets are concatenated into a single MeasurementSet.
    The concatenation is done in an entirely or almost entirely virtual way,
    so hardly any data are copied. It makes the command very fast and hardly
    any extra disk space is needed.

    The MSs can be concatenated in time or frequency (spectral windows).
    If concatenated in time, no indices need to be updated and the
    concatenation is done in a single step.

    If spectral windows are concatenated, tThe data-description-ids and
    spectral-window-ids in the resulting MS and its subtables are updated
    to make them unique.
    The spectral concatenation is done in two steps and results in two MSs:

    1. The input MSs are virtually concatenated resulting in the
       MeasurementSet `<newname>_CONCAT`.
    2. The MeasurementSet <newname> is created. It references all columns
       in `<newname>_CONCAT` with the exception of the DATA_DESC_ID column.
       This column is copied and updated to make the ids correct.
       Furthermore the MS contains a copy of all subtables (with the exception
       of SORTED_TABLE), where the DATA_DESCRIPTION and SPECTRAL_WINDOW
       subtables are the concatenation of those subtables in the input MSs.
       The ids in the resulting subtables are updated.

    The FEED, FREQ_OFFSET, SOURCE, and SYSCAL subtables also have a
    SPECTRAL_WINDOW_ID column. Currently these subtables are not concatenated
    nor are their ids updated.

    `names`
      A sequence containing the names of the MeasurementSets to concatenate.
    `newname`
      The name of the resulting MeasurementSet. A MeasurementSet with this
      name followed by `_CONCAT` will also be created (and must be kept).
    `concatTime`
      False means that the spectral windows ids will be adjusted as explained
      above.

    """
    from pyrap.tables import table,taql
    import numpy as np

    if len(names) == 0:
        raise ValueError('No input MSs given')
    # Concatenation in time is straightforward.
    if concatTime:
        t = table(names[0])
        if 'SYSCAL' in t.fieldnames():
            tn = table(names, concatsubtables='SYSCAL')
        else:
            tn = table(names)
        t.close()
        tn.rename (newname)
        return
    # First concatenate the given tables as another table.
    # The SPECTRAL_WINDOW and DATA_DESCRIPTION subtables are concatenated
    # and changed later.
    # Those subtables cannot be concatenated here, because the deep copy of
    # them fails due to the rename of the main table.
    tn = table(names)
    tdesc = tn.getdesc()
    tn.rename (newname + '_CONCAT')
    tn.flush()
    # Now create a table where all columns forward to the concatenated table,
    # but create a stored column for the data description id, because it has
    # to be changed.
    # The new column is filled at the end.
    tnew = table(newname, tdesc, nrow=tn.nrows(), dminfo={'1':{'TYPE':'ForwardColumnEngine', 'NAME':'ForwardData', 'COLUMNS':tn.colnames(), 'SPEC':{'FORWARDTABLE':tn.name()}}})
    # Remove the DATA_DESC_ID column and recreate it in a stored way.
    tnew.removecols ('DATA_DESC_ID')
    tnew._addcols (pt.tableutil.maketabdesc(pt.tableutil.makecoldesc('DATA_DESC_ID', tdesc['DATA_DESC_ID'])),
                  dminfo={'TYPE':'IncrementalStMan', 'NAME':'DDID', 'SPEC':{}})
    tnew._makerow()

    # Copy the table keywords.
    keywords = tn.getkeywords()
    tnew.putkeywords (keywords)
    # Copy all column keywords.
    for col in tn.colnames():
        tnew.putcolkeywords (col, tn.getcolkeywords(col))
    # Make a deep copy of all subtables (except SORTED_TABLE).
    for key in keywords:
        if key != 'SORTED_TABLE':
            val = keywords[key]
            if isinstance(val, str):
                tsub = table(val, ack=False)
                tsubn = tsub.copy (newname + '/' + key, deep=True)
                tnew.putkeyword (key, tsubn)
    tnew.flush()
    # Now we have to take care that the subbands are numbered correctly.
    # The DATA_DESCRIPTION and SPECTRAL_WINDOW subtables are concatenated.
    # The ddid in the main table and spwid in DD subtable have to be updated.
    tnewdd  = table(tnew.getkeyword('DATA_DESCRIPTION'), readonly=False, ack=False)
    tnewspw = table(tnew.getkeyword('SPECTRAL_WINDOW'), readonly=False, ack=False)
    nrdd   = 0
    nrspw  = 0
    nrmain = 0
    useChanSel = True
    for name in names:
        t = table(name, ack=False)
        tdd  = table(t.getkeyword('DATA_DESCRIPTION'), ack=False)
        tspw = table(t.getkeyword('SPECTRAL_WINDOW'), ack=False)
        # The first table already has its subtable copied.
        # Append the subtables of the other ones.
        if nrdd > 0:
            tnewdd.addrows (tdd.nrows())
            for i in range(tdd.nrows()):
                tnewdd[nrdd+i] = tdd[i]        # copy row i
            tnewspw.addrows (tspw.nrows())
            for i in range(tspw.nrows()):
                tnewspw[nrspw+i] = tspw[i]
        tnewdd.putcol ('SPECTRAL_WINDOW_ID',
                       tdd.getcol('SPECTRAL_WINDOW_ID') + nrspw,
                       nrdd, tdd.nrows())
        tnew.putcol ('DATA_DESC_ID',
                     t.getcol('DATA_DESC_ID') + nrdd,
                     nrmain, t.nrows())
        nrdd += tdd.nrows()
        nrspw += tspw.nrows()
        nrmain += t.nrows()
    # Overwrite keyword CHANNEL_SELECTION.
    if 'MODEL_DATA' in tnew.colnames():
        if 'CHANNEL_SELECTION' in tnew.colkeywordnames('MODEL_DATA'):
            tnew.removecolkeyword ('MODEL_DATA', 'CHANNEL_SELECTION')
            # Define the CHANNEL_SELECTION keyword containing the channels of
            # all spectral windows.
            tspw = table(tnew.getkeyword('SPECTRAL_WINDOW'), ack=False)
            nchans = tspw.getcol('NUM_CHAN')
            chans = [[0,nch] for nch in nchans]
            tnew.putcolkeyword ('MODEL_DATA', 'CHANNEL_SELECTION',
                                np.int32(chans))
    tnew.flush(True)

if __name__ == '__main__':
    descriptiontext = "Perform virtual concatenation.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_files', help='list of ms files to concatenate')
    parser.add_argument('outfile', help='output filename')
    parser.add_argument('-c', '--clobber', help='overwrite existing outfile?', type=bool, default=True)

    args = parser.parse_args()
    main(args.ms_files, args.outfile, clobber=args.clobber)
