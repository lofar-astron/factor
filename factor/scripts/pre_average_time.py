#! /usr/bin/env python
"""
Script to pre-average data using a sliding Gaussian kernel in time
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import glob
import sys
import os
import itertools
import pickle
from scipy.ndimage.filters import gaussian_filter1d as gfilter
import casacore.tables as pt
import lofar.parmdb
from astropy.stats import median_absolute_deviation


def main(ms_input, parmdb_input, input_colname, output_data_colname, output_weights_colname,
    target_rms_rad, baseline_file, minutes_per_block=10.0, verbose=True):
    """
    Pre-average data using a sliding Gaussian kernel in time

    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile
    parmdb_input : list or str
        List of parmdb filenames, or string with list, or path to a mapfile
        The resulting list from parmdb_input must match the one from ms_input
    input_colname : str
        Name of the column in the MS from which the data are read
    output_data_colname : str
        Name of the column in the MS into which the averaged data are written
    output_weights_colname : str
        Name of the column in the MS into which the averaged data weights are written
    target_rms_rad : float (str)
        The target RMS for the phase noise in the input parmDBs. (Or whatever???)
    baseline_file : str
        Filename of pickled baseline lengths

    """

    # convert input to needed types
    ms_list = input2strlist(ms_input)
    parmdb_list = input2strlist(parmdb_input)
    verbose = input2bool(verbose)
    if len(ms_list) != len(parmdb_list):
        raise ValueError('pre_average_time: Length of MS-list ({0}) and length of parmdb-list ({1}) differ.'.format(len(ms_list),len(parmdb_list)))

    if type(target_rms_rad) is str:
        target_rms_rad = float(target_rms_rad)
    if os.path.exists(baseline_file):
        f = open(baseline_file, 'r')
        baseline_dict = pickle.load(f)
        f.close()
    else:
        print('Cannot find baseline_file. Exiting...')
        sys.exit(1)

    # Iterate through time chunks and find the lowest ionfactor
    start_times = []
    end_times = []
    ionfactors = []
    if verbose:
        print('Determining ionfactors...')
    for msind in xrange(len(ms_list)):
        tab = pt.table(ms_list[msind], ack=False)
        start_time = tab[0]['TIME']
        end_time = tab[-1]['TIME']
        remaining_time = end_time - start_time # seconds
        start_times.append(start_time)
        end_times.append(end_time)
        tab.close()
        t_delta = minutes_per_block * 60.0 # seconds
        t1 = 0.0
        while remaining_time > 0.0:
            if remaining_time < 1.5 * t_delta:
                # If remaining time is too short, just include it all in this chunk
                t_delta = remaining_time + 10.0
            remaining_time -= t_delta

            # Find ionfactor for this period
            ionfactors.append(find_ionfactor(parmdb_list[msind], baseline_dict, t1+start_time,
                                             t1+start_time+t_delta, target_rms_rad=target_rms_rad))
            if verbose:
                print('    ionfactor (for timerange {0}-{1} sec) = {2}'.format(t1,
                      t1+t_delta, ionfactors[-1]))
            t1 += t_delta

    sorted_ms_tuples = sorted(zip(start_times,end_times,range(len(ms_list)),ms_list))
    sorted_ms_dict = { 'msnames' :[ms for starttime,endtime,index,ms in sorted_ms_tuples],
                       'starttimes' : [starttime for starttime,endtime,index,ms in sorted_ms_tuples],
                       'endtimes' : [endtime for starttime,endtime,index,ms in sorted_ms_tuples] }

    # Do pre-averaging using lowest ionfactor
    ionfactor_min = min(ionfactors)
    if verbose:
        print('Using ionfactor = {}'.format(ionfactor_min))
        print('Averaging...')
    BLavg_multi(sorted_ms_dict, baseline_dict, input_colname, output_data_colname,
        output_weights_colname, ionfactor_min)


def find_ionfactor(parmdb_file, baseline_dict, t1, t2, target_rms_rad=0.2):
    """
    Finds ionospheric scaling factor
    """
    pdb_in = lofar.parmdb.parmdb(parmdb_file)
    parms = pdb_in.getValuesGrid('*')

    # Filter any stations not in both the instrument table and the ms
    stations_pbd = set([s.split(':')[-1] for s in pdb_in.getNames()])
    stations_ms = set([s for s in baseline_dict.itervalues() if type(s) is str])
    stations = sorted(list(stations_pbd.intersection(stations_ms)))

    # Select long baselines only (BL > 10 km), as they will set the ionfactor scaling
    ant1 = []
    ant2 = []
    dist = []
    min_length = 10.0
    for k, v in baseline_dict.iteritems():
        if type(v) is not str and '-' in k:
            if v > min_length:
                s1 = k.split('-')[0]
                s2 = k.split('-')[1]
                s1_name = baseline_dict[s1]
                s2_name = baseline_dict[s2]
                if s1_name in stations and s2_name in stations:
                    ant1.append(s1_name)
                    ant2.append(s2_name)
                    dist.append(v)

    # Find correlation times
    rmstimes = []
    dists = []
    freq = None
    for a1, a2, d in zip(ant1, ant2, dist):
        if freq is None:
            freq = np.copy(parms['Gain:0:0:Phase:{}'.format(a1)]['freqs'])[0]
            times = np.copy(parms['Gain:0:0:Phase:{}'.format(a1)]['times'])
            time_ind = np.where((times >= t1) & (times < t2))[0]
            timepersolution = np.copy(parms['Gain:0:0:Phase:{}'.format(a1)]['timewidths'])[0]
        ph1 = np.copy(parms['Gain:0:0:Phase:{}'.format(a1)]['values'])[time_ind]
        ph2 = np.copy(parms['Gain:0:0:Phase:{}'.format(a2)]['values'])[time_ind]

        # Filter flagged solutions
        good = np.where((~np.isnan(ph1)) & (~np.isnan(ph2)))[0]
        if len(good) == 0:
            continue

        rmstime = None
        ph = unwrap_fft(ph2[good] - ph1[good])

        step = 1
        for i in range(1, len(ph)/2, step):
            p1 = ph[i:]
            p2 = ph[:-i]
            mean = np.mean(p1-p2)
            if mean > target_rms_rad:
                rmstime = i
                break
        if rmstime is None:
            rmstime = len(ph)/2
        rmstimes.append(rmstime)
        dists.append(d)

    # Find the mean ionfactor assuming that the correlation time goes as
    # t_corr ~ 1/sqrt(BL). The ionfactor is defined in BLavg() as:
    #
    #     ionfactor = (t_corr / 30.0 sec) / ( np.sqrt((25.0 / dist_km)) * (freq_hz / 60.e6) )
    #
    ionfactor = np.mean(np.array(rmstimes) / 30.0 / (np.sqrt(25.0 / np.array(dists))
        * freq / 60.0e6)) * timepersolution

    return ionfactor


def BLavg_multi(sorted_ms_dict, baseline_dict, input_colname, output_data_colname,
        output_weights_colname, ionfactor, clobber=True, maxgap_sec=1800,
        check_files = True):
    """
    Averages data using a sliding Gaussian kernel on the weights
    """

    #### sort msnames into groups with gaps < maxgap_sec
    nfiles = len(sorted_ms_dict['msnames'])
    ms_groups = []
    newgroup = []
    for msindex in xrange(nfiles):
        if msindex+1 == nfiles or sorted_ms_dict['starttimes'][msindex+1] > sorted_ms_dict['endtimes'][msindex] + maxgap_sec:
            newgroup.append(sorted_ms_dict['msnames'][msindex])
            ms_groups.append(newgroup)
            newgroup = []
        else:
            newgroup.append(sorted_ms_dict['msnames'][msindex])

    print "BLavg_multi: Working on",len(ms_groups),"groups of measurement sets."
    #### loop over all groups
    msindex = 0
    for ms_names in ms_groups:
        ### collect data from all files in this group
        freqtab = pt.table(ms_names[0] + '::SPECTRAL_WINDOW', ack=False)
        freq = freqtab.getcell('REF_FREQUENCY',0)
        freqtab.close()
        timepersample = None
        ant1_list        = []
        ant2_list        = []
        all_time_list    = []
        all_data_list    = []
        all_weights_list = []
        all_flags_list   = []
        for msfile in ms_names:
            if not os.path.exists(msfile):
                print("Cannot find MS file: {0}.".format(msfile))
                sys.exit(1)
            # open input/output MS
            ms = pt.table(msfile, readonly=True, ack=False)
            if check_files:
                freqtab = pt.table(msfile + '::SPECTRAL_WINDOW', ack=False)
                if freqtab.getcell('REF_FREQUENCY',0) != freq:
                    print("Different REF_FREQUENCYs: {0} and: {1} in {2}.".format(freq,freqtab.getcell('REF_FREQUENCY',0),msfile))
                    sys.exit(1)
                freqtab.close()
            #wav = 299792458. / freq
            if timepersample is None:
                timepersample = ms.getcell('INTERVAL',0)
            elif check_files:
                if timepersample != ms.getcell('INTERVAL',0):
                    print("Different INTERVALs: {0} and: {1} in {2}.".format(timepersample,ms.getcell('INTERVAL',0),msfile))
                    sys.exit(1)
            all_time_list.append( ms.getcol('TIME_CENTROID') )
            ant1_list.append( ms.getcol('ANTENNA1') )
            ant2_list.append( ms.getcol('ANTENNA2') )
            all_data_list.append( ms.getcol(input_colname) )
            all_weights_list.append( ms.getcol('WEIGHT_SPECTRUM') )
            all_flags_list.append( ms.getcol('FLAG') )

            all_flags_list[-1][ np.isnan(all_data_list[-1]) ] = True # flag NaNs
            all_weights_list[-1] = all_weights_list[-1] * ~all_flags_list[-1] # set weight of flagged data to 0

            # Check that all NaNs are flagged
            if np.count_nonzero(np.isnan(all_data_list[-1][~all_flags_list[-1]])) > 0:
                logging.error('NaNs in unflagged data in {0}!'.format(msfile))
                sys.exit(1)

        ### iteration on baseline combination
        for ant in itertools.product(set(ant1_list[0]), set(ant2_list[0])):
            if ant[0] >= ant[1]:
                continue
            sel_list = []
            weights_list = []
            data_list = []
            # select data from all MSs
            for msindex in xrange(len(ms_names)):
                sel1 = np.where(ant1_list[msindex] == ant[0])[0]
                sel2 = np.where(ant2_list[msindex] == ant[1])[0]
                sel_list.append( sorted(list(frozenset(sel1).intersection(sel2))) )

                # # get weights and data
                # weights_list.append( all_weights[sel_list[msindex],:,:] )
                # data_list.append( all_data[sel_list[msindex],:,:] )
            # combine data and weights into one array
            data = all_data_list[0][sel_list[0],:,:]
            weights = all_weights_list[0][sel_list[0],:,:]
            fillshape = list(data.shape)
            startidx = [0]
            endidx = [data.shape[0]]
            for msindex in xrange(1,len(ms_names)):
                #pad gap between obs
                filltimes = np.arange(np.max(all_time_list[msindex-1]),np.min(all_time_list[msindex]),timepersample)
                fillshape[0] = len(filltimes)
                data = np.concatenate( (data,np.zeros(fillshape)), axis=0 )
                weights = np.concatenate( (weights,np.zeros(fillshape)), axis=0  )
                startidx.append(data.shape[0])
                data = np.concatenate( (data,all_data_list[msindex][sel_list[msindex],:,:]), axis=0  )
                weights = np.concatenate( (weights,all_weights_list[msindex][sel_list[msindex],:,:]), axis=0  )
                endidx.append(data.shape[0])

            # compute the FWHM
            dist = baseline_dict['{0}-{1}'.format(ant[0], ant[1])]
            stddev = 30.0 * ionfactor * np.sqrt((25.0 / dist)) * (freq / 60.e6) # in sec
            stddev = stddev/timepersample # in samples

            # Multiply every element of the data by the weights, convolve both
            # the scaled data and the weights, and then divide the convolved data
            # by the convolved weights (translating flagged data into weight=0).
            # That's basically the equivalent of a running weighted average with
            # a Gaussian window function.

            # weigth data and set bad data to 0 so nans do not propagate
            data = np.nan_to_num(data*weights)

            # smear weighted data and weights
            dataR = gfilter(np.real(data), stddev, axis=0)
            dataI = gfilter(np.imag(data), stddev, axis=0)
            weights = gfilter(weights, stddev, axis=0)

            # re-create data
            data = (dataR + 1j * dataI)
            data[(weights != 0)] /= weights[(weights != 0)] # avoid divbyzero
            for msindex in xrange(len(ms_names)):
                all_data_list[msindex][sel_list[msindex],:,:] = data[startidx[msindex]:endidx[msindex],:,:]
                all_weights_list[msindex][sel_list[msindex],:,:] = weights[startidx[msindex]:endidx[msindex],:,:]

        ### write the data back to the files
        for msindex in xrange(len(ms_names)):
            ms = pt.table(ms_names[msindex], readonly=False, ack=False)
            # Add the output columns if needed
            if output_data_colname not in ms.colnames():
                desc = ms.getcoldesc(input_colname)
                desc['name'] = output_data_colname
                ms.addcols(desc)
            if output_weights_colname not in ms.colnames():
                desc = ms.getcoldesc('WEIGHT_SPECTRUM')
                desc['name'] = output_weights_colname
                ms.addcols(desc)

            ms.putcol(output_data_colname, all_data_list[msindex])
            ms.putcol('FLAG', all_flags_list[msindex]) # this saves flags of nans, which is always good
            ms.putcol(output_weights_colname, all_weights_list[msindex])
            ms.close()
        print "BLavg_multi: Finished one group of measurement sets."


def unwrap_fft(phase, iterations=3):
    """
    Unwrap phase using Fourier techniques.

    For details, see:
    Marvin A. Schofield & Yimei Zhu, Optics Letters, 28, 14 (2003)

    Keyword arguments:
    phase -- array of phase solutions
    iterations -- number of iterations to perform
    """
    puRadius=lambda x : np.roll( np.roll(
          np.add.outer( np.arange(-x.shape[0]/2+1,x.shape[0]/2+1)**2.0,
                        np.arange(-x.shape[1]/2+1,x.shape[1]/2+1)**2.0 ),
          x.shape[1]/2+1,axis=1), x.shape[0]/2+1,axis=0)+1e-9

    idt,dt=np.fft.ifft2,np.fft.fft2
    puOp=lambda x : idt( np.where(puRadius(x)==1e-9,1,puRadius(x)**-1.0)*dt(
          np.cos(x)*idt(puRadius(x)*dt(np.sin(x)))
         -np.sin(x)*idt(puRadius(x)*dt(np.cos(x))) ) )

    def phaseUnwrapper(ip):
       mirrored=np.zeros([x*2 for x in ip.shape])
       mirrored[:ip.shape[0],:ip.shape[1]]=ip
       mirrored[ip.shape[0]:,:ip.shape[1]]=ip[::-1,:]
       mirrored[ip.shape[0]:,ip.shape[1]:]=ip[::-1,::-1]
       mirrored[:ip.shape[0],ip.shape[1]:]=ip[:,::-1]

       return (ip+2*np.pi*
             np.round((puOp(mirrored).real[:ip.shape[0],:ip.shape[1]]-ip)
             /2/np.pi))

    i = 0
    if iterations < 1:
        interations = 1
    while i < iterations:
        i += 1
        phase = phaseUnwrapper(phase)

    return phase


def input2bool(invar):
    if isinstance(invar, bool):
        return invar
    elif isinstance(invar, str):
        if invar.upper() == 'TRUE' or invar == '1':
            return True
        elif invar.upper() == 'FALSE' or invar == '0':
            return False
        else:
            raise ValueError('input2bool: Cannot convert string "'+invar+'" to boolean!')
    elif isinstance(invar, int) or isinstance(invar, float):
        return bool(invar)
    else:
        raise TypeError('input2bool: Unsupported data type:'+str(type(invar)))


def input2strlist(invar):
    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
        else:
            map_in = DataMap.load(invar)
            map_in.iterator = DataMap.SkipIterator
            str_list = []
            for fname in map_in:
                if fname.startswith('[') and fname.endswith(']'):
                    for f in fname.strip('[]').split(','):
                        str_list.append(f.strip(' \'\"'))
                else:
                    str_list.append(fname.strip(' \'\"'))
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return str_list


if __name__ == '__main__':
    descriptiontext = "Pre-average data.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_file_pattern', help='Glob-able filename-pattern of input datasets')
    parser.add_argument('parmdb_file_pattern', help='Glob-able filename-pattern of input direction-independent selfcal instrument parmdbs')
    parser.add_argument('input_colname', help='Name of input column to pre-average')
    parser.add_argument('output_data_colname', help='Name of output column')
    parser.add_argument('output_weights_colname', help='Name of output column')
    parser.add_argument('target_rms', help='Target rms in Jy/beam')
    args = parser.parse_args()

    ms_input = glob.glob(args.ms_file_pattern)
    parmdb_input = glob.glob(args.parmdb_file_pattern)

    main(ms_input, parmdb_input, args.input_colname, args.output_data_colname,
        args.output_weights_colname, args.target_rms)
