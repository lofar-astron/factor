#! /usr/bin/env python
"""
Script to convert selfcal solutions to single gain table
"""
import argparse
from argparse import RawTextHelpFormatter
import lofar.parmdb as lp
import shutil
import numpy as np
import sys
import os


def main(fast_parmdb, slow_parmdb, output_file, freqstep=1, preapply_parmdb=None,
    scratch_dir=None):
    """
    Converts multiple selfcal tables to single gain table

    Parameters
    ----------
    fast_parmdb : str
        File with fast phase (TEC and CommonScalarPhase) solutions
    fast_parmdb : str
        File with slow gain solutions
    output_file : str
        Output filename
    freqstep : int
        Frequency step to divide up frequency width of solutions
    preapply_parmdb : str
        File with combined fast phase (TEC and CommonScalarPhase) and slow phase
        solutions for pre-application
    scratch_dir : str, optional
        Scratch directory for temp storage

    """
    freqstep = int(freqstep)

    # Copy to scratch directory if specified
    if scratch_dir is not None:
        fast_parmdb_orig = fast_parmdb
        fast_parmdb = os.path.join(scratch_dir, os.path.basename(fast_parmdb_orig))
        slow_parmdb_orig = slow_parmdb
        slow_parmdb = os.path.join(scratch_dir, os.path.basename(slow_parmdb_orig))
        output_file_orig = output_file
        output_file = os.path.join(scratch_dir, os.path.basename(output_file_orig))
        shutil.copytree(fast_parmdb_orig, fast_parmdb)
        shutil.copytree(slow_parmdb_orig, slow_parmdb)

    fast_pdb = lp.parmdb(fast_parmdb)
    fast_soldict = fast_pdb.getValuesGrid('*')
    slow_pdb = lp.parmdb(slow_parmdb)
    slow_soldict = slow_pdb.getValuesGrid('*')
    output_pdb = lp.parmdb(output_file, create=True)
    if preapply_parmdb is not None:
        preapply_pdb = lp.parmdb(preapply_parmdb)
        preapply_soldict = preapply_pdb.getValuesGrid('*')

    # Get various quantities over which we must iterate
    station_names = list(set([s.split(':')[-1] for s in fast_pdb.getNames()]))
    fast_times = fast_soldict['CommonScalarPhase:{s}'.format(s=station_names[0])]['times']
    fast_timewidths = fast_soldict['CommonScalarPhase:{s}'.format(s=station_names[0])]['timewidths']
    fast_timestep = np.mean(fast_timewidths)

    if preapply_parmdb is not None:
        fast_times_preapply = preapply_soldict['Gain:0:0:Phase:{s}'.format(s=station_names[0])]['times']
        fast_timewidths_preapply = preapply_soldict['Gain:0:0:Phase:{s}'.format(s=station_names[0])]['timewidths']
        fast_timestep_preapply = np.mean(fast_timewidths_preapply)

    slow_freqs = slow_soldict['Gain:0:0:Real:{s}'.format(s=station_names[0])]['freqs']
    slow_freqwidths = slow_soldict['Gain:0:0:Real:{s}'.format(s=station_names[0])]['freqwidths']
    slow_freqstep = np.mean(slow_freqwidths)
    if preapply_parmdb is not None:
        slow_freqs_preapply = preapply_soldict['Gain:0:0:Phase:{s}'.format(s=station_names[0])]['freqs']
        slow_freqwidths_preapply = preapply_soldict['Gain:0:0:Phase:{s}'.format(s=station_names[0])]['freqwidths']
        slow_freqstep_preapply = np.mean(slow_freqwidths_preapply)

    key_names = fast_soldict.keys()
    if any('Gain:0:1:' in s for s in key_names):
        pol_list = ['0:0', '1:1', '0:1', '1:0']
    else:
        pol_list = ['0:0', '1:1']

    # Check preapply to make sure we use the finest grid
    if preapply_parmdb is not None:
        if slow_freqstep_preapply < slow_freqstep:
            slow_freqs = slow_freqs_preapply
            slow_freqwidths = slow_freqwidths_preapply
        if fast_timestep_preapply < fast_timestep:
            fast_times = fast_times_preapply
            fast_timewidths = fast_timewidths_preapply

    # Determine final time and frequency grid. This step is not needed if a
    # preapply_parmdb is specified, as it will already be made on the final grids
    if freqstep > 1 and preapply_parmdb is None:
        final_freqs = []
        final_freqwidths = []
        for freq, freqwidth in zip(slow_freqs, slow_freqwidths):
            final_freqwidth = freqwidth / freqstep
            low_freq = freq - freqwidth / 2
            for i in range(freqstep):
                final_freqs.append(low_freq + final_freqwidth * (i + 0.5))
                final_freqwidths.append(final_freqwidth)
        final_freqs = np.array(final_freqs)
        final_freqwidths = np.array(final_freqwidths)
    else:
        final_freqs = slow_freqs
        final_freqwidths = slow_freqwidths

    # Identify any gaps in time (frequency gaps are not allowed), as we need to handle
    # each section separately if gaps are present
    delta_times = fast_times[1:] - fast_times[:-1]
    gaps = np.where(delta_times > fast_timewidths[:-1]*2.)
    gaps_ind = gaps[0] + 1
    gaps_ind = np.append(gaps_ind, np.array([len(fast_times)]))

    # Get values on the final time and frequency grid
    g_start = 0
    for g in gaps_ind:
        if preapply_parmdb is not None:
            preapply_soldict = preapply_pdb.getValues('*', final_freqs, final_freqwidths,
                fast_times[g_start:g], fast_timewidths[g_start:g], asStartEnd=False)
        fast_soldict = fast_pdb.getValues('*', final_freqs, final_freqwidths,
            fast_times[g_start:g], fast_timewidths[g_start:g], asStartEnd=False)
        slow_soldict = slow_pdb.getValues('*', final_freqs, final_freqwidths,
            fast_times[g_start:g], fast_timewidths[g_start:g], asStartEnd=False)

        # Add various phase and amp corrections together
        for station in station_names:
            try:
                fast_phase = np.copy(fast_soldict['CommonScalarPhase:{s}'.format(s=station)]['values'])
                tec = np.copy(fast_soldict['TEC:{s}'.format(s=station)]['values'])
                tec_phase =  -8.44797245e9 * tec / final_freqs

                for pol in pol_list:
                    slow_real = np.copy(slow_soldict['Gain:'+pol+':Real:{s}'.format(s=station)]['values'])
                    slow_imag = np.copy(slow_soldict['Gain:'+pol+':Imag:{s}'.format(s=station)]['values'])
                    slow_amp = np.sqrt((slow_real**2) + (slow_imag**2))
                    slow_phase = np.arctan2(slow_imag, slow_real)

                    total_amp = slow_amp
                    if preapply_parmdb is not None:
                        fast_phase_preapply = np.copy(preapply_soldict['Gain:'+pol+':Phase:{s}'.format(s=station)]['values'])
                        total_phase = np.mod(fast_phase + tec_phase + slow_phase +
                            fast_phase_preapply + np.pi, 2*np.pi) - np.pi

                        # Identify zero phase solutions and set the corresponding entries in total_phase and total_amp to NaN
                        total_phase = np.where(fast_phase_preapply == 0.0, np.nan, total_phase)
                        total_amp = np.where(fast_phase_preapply == 0.0, np.nan, total_amp)
                    else:
                        total_phase = np.mod(fast_phase + tec_phase + slow_phase + np.pi, 2*np.pi) - np.pi

                    # Identify zero phase solutions and set the corresponding entries in total_phase and total_amp to NaN
                    total_phase = np.where(np.logical_or(fast_phase == 0.0, tec_phase == 0.0), np.nan, total_phase)
                    total_amp = np.where(np.logical_or(fast_phase == 0.0, tec_phase == 0.0), np.nan, total_amp)

                    output_pdb.addValues('Gain:'+pol+':Phase:{}'.format(station),
                        total_phase, final_freqs, final_freqwidths,
                        fast_times[g_start:g], fast_timewidths[g_start:g], asStartEnd=False)
                    output_pdb.addValues('Gain:'+pol+':Ampl:{}'.format(station),
                        total_amp, final_freqs, final_freqwidths,
                        fast_times[g_start:g], fast_timewidths[g_start:g], asStartEnd=False)
            except KeyError:
                # The station in question does not have solutions for this
                # time period
                continue
        g_start = g

    # Write values
    output_pdb.flush()

    # Copy output to original path and delete copies if scratch directory is specified
    if scratch_dir is not None:
        if os.path.exists(output_file_orig):
            shutil.rmtree(output_file_orig)
        shutil.copytree(output_file, output_file_orig)
        shutil.rmtree(output_file)
        shutil.rmtree(fast_parmdb)
        shutil.rmtree(slow_parmdb)


if __name__ == '__main__':
    descriptiontext = "Converts multiple selfcal tables to single gain table.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fast_selfcal_parmdb', help='name of the parmdb with fast solutions')
    parser.add_argument('slow_selfcal_parmdb', help='name of the parmdb with slow solutions')
    parser.add_argument('output_file', help='name of the output file')
    parser.add_argument('--preapply_parmdb',type=str,default=None,help='preapply parmDB from another facet')
    args = parser.parse_args()

    main(args.fast_selfcal_parmdb, args.slow_selfcal_parmdb, args.output_file, preapply_parmdb=args.preapply_parmdb)
