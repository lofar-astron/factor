#! /usr/bin/env python
"""
Script to check a list of MS files for missing frequencies
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import sys


def main(ms_list):
    """
    Check a list of MS files for missing frequencies

    Parameters
    ----------
    ms_list : list
        List of MS filenames, in order of increasing frequency

    Returns
    -------
    result : dict
        Dict with list of MS filenames with dummy filenames inserted in the gaps

    """
    if type(ms_list) is str:
        ms_list = [f.strip() for f in ms_list.strip('[]').split(',')]

    freqs = []
    for i, ms in enumerate(ms_list):
        # Get the frequency info
        sw = pt.table(ms+'::SPECTRAL_WINDOW', ack=False)
        if i == 0:
            freq_width = sw.col('TOTAL_BANDWIDTH')[0]
        freqs.append(sw.col('REF_FREQUENCY')[0])
        sw.close()

    # Find gaps, if any
    missing_bands = []
    for i, (freq1, freq2) in enumerate(zip(freqs[:-1], freqs[1:])):
        ngap = int(round((freq2 - freq1)/freq_width))
        missing_bands.extend([i + j + 1 for j in range(ngap)])

    for m in reversed(missing_bands):
        ms_list.insert(m, 'dummy.ms')

    return {'padded_list': '[{}]'.format(','.join(ms_list))}


if __name__ == '__main__':
    descriptiontext = "Add/subtract columns (column_out = column1 +/- column2).\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_list', help='list of MS files')
    args = parser.parse_args()

    main(args.ms_list)
