#! /usr/bin/env python
"""
Script to check a list of MS files for missing frequencies
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import numpy as np
import sys


def main(ms_list, num_groups=1):
    """
    Check a list of MS files for missing frequencies

    Parameters
    ----------
    ms_list : list
        List of MS filenames, in order of increasing frequency
    max_per_group : int
        Maximum number of bands per group

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
        missing_bands.extend([i + j + 1 for j in range(ngap-1)])

    # Insert dummy entries in gaps
    for m in reversed(missing_bands):
        ms_list.insert(m, 'dummy.ms')

    # Group and convert lists to strings. If there are more than one group,
    # the returned mapfile will have the format:
    # [[ms1, ms2, ms3]||[ms4, ms5, ms6]]
    ms_groups = np.array_split(ms_list, num_groups)
    padded_list = []
    for ms_group in ms_groups:
        padded_list.append('[{}]'.format(','.join(ms_group)))
    padded_list_joined = '||'.join(padded_list)

    return {'padded_list': padded_list}


if __name__ == '__main__':
    descriptiontext = "Check a list of MS files for missing frequencies.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_list', help='list of MS files')
    args = parser.parse_args()

    main(args.ms_list)
