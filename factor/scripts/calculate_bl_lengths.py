#! /usr/bin/env python
"""
Calculates the baseline lengths in km for a MS
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import sys
import os
import itertools
import casacore.tables as pt
import pickle


def main(ms_input, output_file):
    """
    Calculates the baseline lengths in km for a MS

    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile. The MS
        files should all have the same frequency
    output_file : str
        Filename of output pickle file

    """
    ms_list = input2strlist(ms_input)

    print('Calculating baseline lengths...')
    baseline_dict = get_baseline_lengths(ms_list[0])

    with open(output_file, 'wb') as f:
        pickle.dump(baseline_dict, f)


def input2strlist(invar):
    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return str_list


def get_baseline_lengths(ms):
    """
    Returns dict of baseline lengths in km for all baselines in input dataset
    """
    anttab = pt.table(ms+'::ANTENNA', ack=False)
    antnames = anttab.getcol('NAME')
    anttab.close()

    t = pt.table(ms, ack=False)
    ant1 = t.getcol('ANTENNA1')
    ant2 = t.getcol('ANTENNA2')
    all_uvw = t.getcol('UVW')
    t.close()

    baseline_dict = {}
    for ant in itertools.product(set(ant1), set(ant2)):
        if ant[0] >= ant[1]:
            continue
        sel1 = np.where(ant1 == ant[0])[0]
        sel2 = np.where(ant2 == ant[1])[0]
        sel = sorted(list(frozenset(sel1).intersection(sel2)))
        uvw = all_uvw[sel, :]
        uvw_dist = np.sqrt(uvw[:, 0]**2 + uvw[:, 1]**2 + uvw[:, 2]**2)
        baseline_dict['{0}'.format(ant[0])] = antnames[ant[0]]
        baseline_dict['{0}'.format(ant[1])] = antnames[ant[1]]
        baseline_dict['{0}-{1}'.format(ant[0], ant[1])] = np.mean(uvw_dist) / 1.e3

    return baseline_dict
