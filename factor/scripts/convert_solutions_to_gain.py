#! /usr/bin/env python
"""
Script to convert selfcal solutions to single gain table
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import lofar.parmdb as lp
import numpy as np
import sys
import os


def main(merged_selfcal_parmdb, output_file, preapply_parmdb=None):
    """
    Converts multiple selfcal tables to single gain table

    Parameters
    ----------
    merged_selfcal_parmdb : str
        File with fast phase (TEC and CommonScalarPhase) and slow gain solutions
    output_file : str
        Output filename
    preapply_parmdb : str
        File with fast phase (TEC and CommonScalarPhase) and slow phase solutions
        for pre-application
    """
    parmdb = lp.parmdb(merged_selfcal_parmdb)
    soldict = parmdb.getValuesGrid('*')
    pdb_out = lp.parmdb(output_file, create=True)
    if preapply_parmdb is not None:
        pdb_preapply = lp.parmdb(preapply_parmdb)
        soldict_preapply = parmdb.getValuesGrid('*')

    # Get various quantities over which we must iterate
    station_names = list(set([s.split(':')[-1] for s in parmdb.getNames()]))
    fast_times = soldict['CommonScalarPhase:{s}'.format(s=station_names[0])]['times']
    fast_timewidths = soldict['CommonScalarPhase:{s}'.format(s=station_names[0])]['timewidths']
    fast_timestep = np.mean(fast_timewidths)

    if preapply_parmdb is not None:
        fast_times_preapply = soldict_preapply['CommonScalarPhase:{s}'.format(s=station_names[0])]['times']
        fast_timewidths_preapply = soldict_preapply['CommonScalarPhase:{s}'.format(s=station_names[0])]['timewidths']
        fast_timestep_preapply = np.mean(fast_timewidths_preapply)

    slow_freqs = soldict['Gain:0:0:Real:{s}'.format(s=station_names[0])]['freqs']
    slow_freqwidths = soldict['Gain:0:0:Real:{s}'.format(s=station_names[0])]['freqwidths']
    slow_freqstep = np.mean(slow_freqwidths)
    if preapply_parmdb is not None:
        slow_freqs_preapply = soldict_preapply['Gain:0:0:Real:{s}'.format(s=station_names[0])]['freqs']
        slow_freqwidths_preapply = soldict_preapply['Gain:0:0:Real:{s}'.format(s=station_names[0])]['freqwidths']
        slow_freqstep_preapply = np.mean(slow_freqwidths_preapply)

    key_names = soldict.keys()
    if any('Gain:0:1:' in s for s in key_names):
        pol_list = ['0:0', '1:1', '0:1', '1:0']
    else:
        pol_list = ['0:0', '1:1']

    # Get values, assuming a fast time grid and a slow freq grid
    if preapply_parmdb is not None:
        # Make sure we use the finest grid
        if slow_freqstep_preapply < slow_freqstep:
            slow_freqs = slow_freqs_preapply
            slow_freqwidths = slow_freqwidths_preapply
        if fast_timestep_preapply < fast_timestep:
            fast_times = fast_times_preapply
            fast_timewidths = fast_timewidths_preapply
        soldict_preapply = pdb_preapply.getValues('*', slow_freqs, slow_freqwidths,
            fast_times, fast_timewidths, asStartEnd=False)
    soldict = parmdb.getValues('*', slow_freqs, slow_freqwidths, fast_times,
        fast_timewidths, asStartEnd=False)

    # Add various phase and amp corrections together
    for station in station_names:
        fast_phase = np.copy(soldict['CommonScalarPhase:{s}'.format(s=station)]['values'])
        tec = np.copy(soldict['TEC:{s}'.format(s=station)]['values'])
        tec_phase =  -8.44797245e9 * tec / slow_freqs

        slow_real_00 = np.copy(soldict['Gain:0:0:Real:{s}'.format(s=station)]['values'])
        slow_imag_00 = np.copy(soldict['Gain:0:0:Imag:{s}'.format(s=station)]['values'])
        slow_amp_00 = np.sqrt((slow_real_00**2) + (slow_imag_00**2))
        slow_phase_00 = np.arctan2(slow_imag_00, slow_real_00)
        slow_real_11 = np.copy(soldict['Gain:1:1:Real:{s}'.format(s=station)]['values'])
        slow_imag_11 = np.copy(soldict['Gain:1:1:Imag:{s}'.format(s=station)]['values'])
        slow_amp_11 = np.sqrt((slow_real_11**2) + (slow_imag_11**2))
        slow_phase_11 = np.arctan2(slow_imag_11, slow_real_11)

        if preapply_parmdb is not None:
            fast_phase_preapply = np.copy(soldict_preapply['CommonScalarPhase:{s}'.format(s=station)]['values'])
            tec_preapply = np.copy(soldict_preapply['TEC:{s}'.format(s=station)]['values'])
            tec_phase_preapply =  -8.44797245e9 * tec_preapply / slow_freqs

            slow_real_00_preapply = np.copy(soldict_preapply['Gain:0:0:Real:{s}'.format(s=station)]['values'])
            slow_imag_00_preapply = np.copy(soldict_preapply['Gain:0:0:Imag:{s}'.format(s=station)]['values'])
            slow_phase_00_preapply = np.arctan2(slow_imag_00_preapply, slow_real_00_preapply)
            slow_real_11_preapply = np.copy(soldict_preapply['Gain:1:1:Real:{s}'.format(s=station)]['values'])
            slow_imag_11_preapply = np.copy(soldict_preapply['Gain:1:1:Imag:{s}'.format(s=station)]['values'])
            slow_phase_11_preapply = np.arctan2(slow_imag_11_preapply, slow_real_11_preapply)

            total_phase_00 = np.mod(fast_phase + tec_phase + slow_phase_00 +
                fast_phase_preapply + tec_phase_preapply + slow_phase_00_preapply + np.pi, 2*np.pi) - np.pi
            total_amp_00 = slow_amp_00
            total_phase_11 = np.mod(fast_phase + tec_phase + slow_phase_11 +
                fast_phase_preapply + tec_phase_preapply + slow_phase_11_preapply + np.pi, 2*np.pi) - np.pi
            total_amp_11 = slow_amp_11
        else:
            total_phase_00 = np.mod(fast_phase + tec_phase + slow_phase_00 + np.pi, 2*np.pi) - np.pi
            total_amp_00 = slow_amp_00
            total_phase_11 = np.mod(fast_phase + tec_phase + slow_phase_11 + np.pi, 2*np.pi) - np.pi
            total_amp_11 = slow_amp_11

        pdb_out.addValues({'Gain:0:0:Phase:{}'.format(station): {'freqs': slow_freqs, 'freqwidths':
            slow_freqwidths, 'times': fast_times, 'timewidths': fast_timewidths, 'values': total_phase_00}})
        pdb_out.addValues({'Gain:0:0:Ampl:{}'.format(station): {'freqs': slow_freqs, 'freqwidths':
            slow_freqwidths, 'times': fast_times, 'timewidths': fast_timewidths, 'values': total_amp_00}})
        pdb_out.addValues({'Gain:1:1:Phase:{}'.format(station): {'freqs': slow_freqs, 'freqwidths':
            slow_freqwidths, 'times': fast_times, 'timewidths': fast_timewidths, 'values': total_phase_11}})
        pdb_out.addValues({'Gain:1:1:Ampl:{}'.format(station): {'freqs': slow_freqs, 'freqwidths':
            slow_freqwidths, 'times': fast_times, 'timewidths': fast_timewidths, 'values': total_amp_11}})

    # Write values
    pdb_out.flush()


if __name__ == '__main__':
    descriptiontext = "Converts multiple selfcal tables to single gain table.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('merged_selfcal_parmdb', help='name of the merged parmdb')
    parser.add_argument('output_file', help='name of the output file')
    args = parser.parse_args()

    main(args.merged_selfcal_parmdb, args.output_file)
