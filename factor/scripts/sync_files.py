#! /usr/bin/env python
"""
Script to sync files
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import numpy as np
import sys
import os


def main(file_from, file_to, use_compression=False):
    """
    Sync a file

    Parameters
    ----------
    file_from : str
        Name of file to from
    file_to : str
        Name of file to copy to
    use_compression : bool, optional
        If True, use Dysco compression to compress copy. NOTE: this options only
        works on DATA, BLAVG_DATA, WEIGHT_SPECTRUM, and BLAVG_WEIGHT_SPECTRUM
        columns and assumes there is no CORRECTED_DATA already

    """
    if file_from.endswith('/'):
        file_from = file_from.strip('/')
    if file_to.endswith('/'):
        file_to = file_to.strip('/')

    destination_dir = os.path.dirname(file_to)
    if len(destination_dir) > 0:
        os.system('/bin/mkdir -p {0}'.format(destination_dir))

    try:
        # Copy files
        if use_compression:
            file_to_tmp = '{}_tmp'.format(file_to)
            if os.path.exists(file_to_tmp):
                os.system('/bin/rm -r {}'.format(file_to_tmp))
            os.system('/bin/cp -r {0} {1}'.format(file_from, file_to_tmp))
        else:
            os.system('/bin/cp -r {0} {1}'.format(file_from, file_to))

        if use_compression:
            # Create compressed version
            t1 = pt.table(file_to_tmp, readonly=False, ack=False)

            # Make compressed versions of DATA, BLAVG_DATA, WEIGHT_SPECTRUM,
            # and BLAVG_WEIGHT_SPECTRUM columns
            for colname_final in ['DATA', 'BLAVG_DATA', 'WEIGHT_SPECTRUM', 'BLAVG_WEIGHT_SPECTRUM']:
                if colname_final in t1.colnames():
                    if colname_final == 'DATA':
                        colname_temp = 'DATA_TEMP'
                        colname_to_copy_from = 'DATA_TEMP'
                        colname_intermediate = 'DATA'
                    elif colname_final == 'BLAVG_DATA':
                        colname_temp = 'DATA_TEMP'
                        colname_to_copy_from = colname_final
                        colname_intermediate = 'DATA'
                    elif colname_final == 'WEIGHT_SPECTRUM':
                        colname_temp = 'WEIGHT_SPECTRUM_TEMP'
                        colname_to_copy_from = 'WEIGHT_SPECTRUM_TEMP'
                        colname_intermediate = 'WEIGHT_SPECTRUM'
                    elif colname_final == 'BLAVG_WEIGHT_SPECTRUM':
                        colname_temp = 'WEIGHT_SPECTRUM_TEMP'
                        colname_to_copy_from = colname_final
                        colname_intermediate = 'WEIGHT_SPECTRUM'
                    t1.renamecol(colname_intermediate, colname_temp)
                    desc = t1.getcoldesc(colname_to_copy_from)
                    desc['name'] = colname_intermediate

                    # Set DyscoStMan to be storage manager for CORRECTED_DATA and WEIGHT_SPECTRUM
                    # We use a visibility bit rate of 16 and truncation of 1.5 sigma to keep the
                    # compression noise below ~ 0.01 mJy, as estimated from Fig 4 of
                    # Offringa (2016). For the weights, we use a bit rate of 12, as
                    # recommended in Sec 4.4 of Offringa (2016)
                    dmi = {
                        'SPEC': {
                            'dataBitCount': np.uint32(16),
                            'distribution': 'Gaussian',
                            'distributionTruncation': 1.5,
                            'normalization': 'RF',
                            'weightBitCount': np.uint32(12)},
                        'NAME': '{}_dm'.format(colname_final),
                        'SEQNR': 1,
                        'TYPE': 'DyscoStMan'}
                    desc['option'] = 1 # make a Direct column
                    t1.addcols(desc, dmi)

                    # Copy data and replace flagged values with NaNs before doing compression
                    data1 = t1.getcol(colname_to_copy_from)
                    flags = t1.getcol('FLAG')
                    flagged = np.where(flags)
                    data1[flagged] = np.NaN
                    t1.putcol(colname_intermediate, data1)
                    if colname_final != colname_intermediate:
                        t1.renamecol(colname_intermediate, colname_final)
                        t1.renamecol(colname_temp, colname_intermediate)
                    else:
                        t1.removecols([colname_temp])

            t1.copy(file_to, deep=True) # write compressed version to output file
            t1.close()
            os.system('rm -r {}'.format(file_to_tmp)) # delete tmp version
    except:
        print('ERROR: could not sync file to {}, likely due to lack of space. '
            'Please make more space available or select a different directory.'.format(destination_dir))
        if os.path.exists(file_to_tmp):
            os.system('rm -r {}'.format(file_to_tmp))
        if os.path.exists(file_to):
            os.system('rm -rf {}'.format(file_to))
        sys.exit(1)


if __name__ == '__main__':
    descriptiontext = "Sync a file.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('file_from', help='name of the image to copy')
    parser.add_argument('file_to', help='loop counter')
    args = parser.parse_args()

    main(args.file_from, args.file_to)





