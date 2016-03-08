"""
Definition of the band class and a few related functions
"""
import os
import sys
import shutil
import logging
import pyrap.tables as pt
import lofar.parmdb
import numpy as np
import multiprocessing
import itertools


class Band(object):
    """
    The Band object contains parameters needed for each band

    Parameters
    ----------
    MSfile : str
        Filename of MS
    factor_working_dir : str
        Full path of working directory
    dirindparmdb : str
        Name of direction-independent instrument parmdb (relative to MSfile)
    skymodel_dirindep : str
        Full path of direction-independent sky model
    local_dir : str
        Path to local scratch directory for temp output. The file is then
        copied to the original output directory
    test_run : bool, optional
        If True, use test image sizes

    """
    def __init__(self, MSfiles, factor_working_dir, dirindparmdb,
        skymodel_dirindep=None, local_dir=None, test_run=False):

        self.files = MSfiles
        self.msnames = [ MS.split('/')[-1] for MS in self.files ]
        self.working_dir = factor_working_dir
        self.dirindparmdbs = [ os.path.join(MS, dirindparmdb) for MS in self.files ]
        self.skymodel_dirindep = skymodel_dirindep
        self.numMS = len(self.files)

        # Get the frequency info and set name
        sw = pt.table(self.files[0]+'::SPECTRAL_WINDOW', ack=False)
        self.freq = sw.col('REF_FREQUENCY')[0]
        self.nchan = sw.col('NUM_CHAN')[0]
        self.chan_freqs_hz = sw.col('CHAN_FREQ')[0]
        self.chan_width_hz = sw.col('CHAN_WIDTH')[0][0]
        sw.close()
        self.name = 'Band_{0:.2f}MHz'.format(self.freq/1e6)
        self.log = logging.getLogger('factor:{}'.format(self.name))
        self.log.debug('Band name is {}'.format(self.name))
        self.chunks_dir = os.path.join(factor_working_dir, 'chunks', self.name)

        # Do some checks
        self.check_freqs()
        self.check_parmdb()

        # Get the field RA and Dec
        obs = pt.table(self.files[0]+'::FIELD', ack=False)
        self.ra = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][0]))
        if self.ra < 0.:
            self.ra = 360.0 + (self.ra)
        self.dec = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][1]))
        obs.close()

        # Get the station diameter
        ant = pt.table(self.files[0]+'::ANTENNA', ack=False)
        self.diam = float(ant.col('DISH_DIAMETER')[0])
        ant.close()

        # Find mean elevation and FOV
        for MS_id in xrange(self.numMS):
            # Add (virtual) elevation column to MS
            try:
                pt.addDerivedMSCal(self.files[MS_id])
            except RuntimeError:
                # RuntimeError indicates column already exists
                pass

            # Calculate mean elevation
            tab = pt.table(self.files[MS_id], ack=False)
            if MS_id == 0:
                global_el_values = tab.getcol('AZEL1', rowincr=10000)[:, 1]
            else:
                global_el_values = np.hstack( (global_el_values, tab.getcol('AZEL1', rowincr=10000)[:, 1]) )
            tab.close()

            # Remove (virtual) elevation column from MS
            pt.removeDerivedMSCal(self.files[MS_id])
        self.mean_el_rad = np.mean(global_el_values)
        sec_el = 1.0 / np.sin(self.mean_el_rad)
        self.fwhm_deg = 1.1 * ((3.0e8 / self.freq) / self.diam) * 180. / np.pi * sec_el

        # Check for SUBTRACTED_DATA_ALL column in original datasets
        self.has_sub_data = True
        self.has_sub_data_new = False
        for MSid in xrange(self.numMS):
            tab = pt.table(self.files[MSid], ack=False)
            if not 'SUBTRACTED_DATA_ALL' in tab.colnames():
                self.log.error('SUBTRACTED_DATA_ALL column not found in file '
                    '{}'.format(self.files[MSid]))
                self.has_sub_data = False
            tab.close()
        if not self.has_sub_data:
            self.log.info('Exiting...')
            sys.exit(1)

        # cut input files into chunks if needed
        chunksize = 2400. # in seconds -> 40min
        self.chunk_input_files(chunksize, dirindparmdb, local_dir=local_dir,
            clobber=True, test_run=test_run)

        # Calculate times and number of samples
        self.sumsamples = 0
        self.minSamplesPerFile = 4294967295  # If LOFAR lasts that many seconds then I buy you a beer.
        self.starttime = np.finfo('d').max
        self.endtime = 0.
        for MSid in xrange(self.numMS):
            tab = pt.table(self.files[MSid], ack=False)
            self.starttime = min(self.starttime,np.min(tab.getcol('TIME')))
            self.endtime = max(self.endtime,np.min(tab.getcol('TIME')))
            for t2 in tab.iter(["ANTENNA1","ANTENNA2"]):
                if (t2.getcell('ANTENNA1',0)) < (t2.getcell('ANTENNA2',0)):
                    self.timepersample = t2.col('TIME')[1] - t2.col('TIME')[0]
                    numsamples = t2.nrows()
                    self.sumsamples += numsamples
                    self.minSamplesPerFile = min(self.minSamplesPerFile,numsamples)
                    break
            tab.close()

        self.log.debug("Using {0} files.".format(len(self.files)))


    def check_parmdb(self):
        """
        Checks the dir-indep instrument parmdb for various problems
        """
        for pdb_id in xrange(self.numMS):
            # Check for special BBS table name "instrument"
            if os.path.basename(self.dirindparmdbs[pdb_id]) == 'instrument':
                self.dirindparmdbs[pdb_id] += '_dirindep'
                if not os.path.exists(self.dirindparmdbs[pdb_id]):
                    if not os.path.exists(os.path.join(self.files[pdb_id], 'instrument')):
                        self.log.critical('Direction-independent instument parmdb not found '
                            'for band {0}'.format(self.files[pdb_id]))
                        sys.exit(1)
                    self.log.warn('Direction-independent instument parmdb for band {0} is '
                        'named "instrument". Copying to "instrument_dirindep" so that BBS '
                        'will not overwrite this table...'.format(self.files[pdb_id]))
                    os.system('cp -r {0} {1}'.format(os.path.join(self.files[pdb_id],
                        'instrument'), self.dirindparmdbs[pdb_id]))
            if not os.path.exists(self.dirindparmdbs[pdb_id]):
                self.log.critical('Direction-independent instrument parmdb "{0}" not found '
                    'for band {1}'.format(self.dirindparmdbs[pdb_id], self.files[pdb_id]))
                sys.exit(1)

            # Check whether there are ampl/phase or real/imag
            try:
                pdb = lofar.parmdb.parmdb(self.dirindparmdbs[pdb_id])
                solname = pdb.getNames()[0]
            except IndexError:
                self.log.critical('Direction-independent instument parmdb appears to be empty '
                            'for band {0}'.format(self.files[pdb_id]))
                sys.exit(1)
            if solname[0:4] != 'Gain':
                self.log.critical('Direction-independent instument parmdb contains not-handled value {0} '
                                  'for band {1}'.format(solname,self.files[pdb_id]))
                sys.exit(1)

            if 'Real' in solname or 'Imag' in solname:
                # Convert real/imag to phasors
                self.log.warn('Direction-independent instument parmdb for band {0} contains '
                    'real/imaginary values. Converting to phase/amplitude...'.format(self.files[pdb_id]))
                self.convert_parmdb_to_phasors_id(pdb_id)
            pdb = False

            # Check that there aren't extra default values in the parmdb, as this
            # confuses DPPP
            pdb = lofar.parmdb.parmdb(self.dirindparmdbs[pdb_id])
            solname = pdb.getNames()[0]
            defvals = pdb.getDefValues()
            for v in defvals:
                if 'Ampl' not in v and 'Phase' not in v:
                    pdb.deleteDefValues(v)
            pdb.flush()


    def convert_parmdb_to_phasors_id(self, pdb_id=0):
        """
        Converts a single instrument parmdb from real/imag to phasors

        Parameters
        ----------
        pdb_id : int
            index of the instrument parmdb to convert
        """
        phasors_parmdb_file = self.dirindparmdbs[pdb_id] + '_phasors'
        pdb_in = lofar.parmdb.parmdb(self.dirindparmdbs[pdb_id])
        pdb_out = lofar.parmdb.parmdb(phasors_parmdb_file, create=True)

        # Check parmdb for non-handled values
        solnames = pdb_in.getNames()
        for name in solnames:
            if name[0:9] != 'Gain:0:0:' and name[0:9] != 'Gain:1:1:':
                self.log.critical('Direction-independent instument parmdb contains not-handled value {0} '
                                  'for band {1}'.format(name,self.files[pdb_id]))
                sys.exit(1)

        # Get station names
        stations = set([s.split(':')[-1] for s in pdb_in.getNames()])

        # Calculate and store phase and amp values for each station
        parms = pdb_in.getValuesGrid('*')
        for i, s in enumerate(stations):
            if i == 0:
                freqs = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['freqs'])
                freqwidths = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['freqwidths'])
                times = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['times'])
                timewidths = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['timewidths'])

            valIm_00 = np.copy(parms['Gain:0:0:Imag:{}'.format(s)]['values'][:, 0])
            valIm_11 = np.copy(parms['Gain:1:1:Imag:{}'.format(s)]['values'][:, 0])
            valRe_00 = np.copy(parms['Gain:0:0:Real:{}'.format(s)]['values'][:, 0])
            valRe_11 = np.copy(parms['Gain:1:1:Real:{}'.format(s)]['values'][:, 0])

            valAmp_00 = np.sqrt((valRe_00**2) + (valIm_00**2))
            valAmp_11 = np.sqrt((valRe_11**2) + (valIm_11**2))
            valPh_00 = np.arctan2(valIm_00, valRe_00)
            valPh_11 = np.arctan2(valIm_11, valRe_11)

            pdb_out.addValues({'Gain:0:0:Phase:{}'.format(s): {'freqs': freqs, 'freqwidths':
                freqwidths, 'times': times, 'timewidths': timewidths, 'values': valPh_00[:,np.newaxis]}})
            pdb_out.addValues({'Gain:1:1:Phase:{}'.format(s): {'freqs': freqs, 'freqwidths':
                freqwidths, 'times': times, 'timewidths': timewidths, 'values': valPh_11[:,np.newaxis]}})
            pdb_out.addValues({'Gain:0:0:Ampl:{}'.format(s): {'freqs': freqs, 'freqwidths':
                freqwidths, 'times': times, 'timewidths': timewidths, 'values': valAmp_00[:,np.newaxis]}})
            pdb_out.addValues({'Gain:1:1:Ampl:{}'.format(s): {'freqs': freqs, 'freqwidths':
                freqwidths, 'times': times, 'timewidths': timewidths, 'values': valAmp_11[:,np.newaxis]}})

        # Write values
        pdb_out.flush()
        pdb_in = False
        pdb_out = False
        self.dirindparmdbs[pdb_id] = phasors_parmdb_file


    def check_freqs(self):
        """
        Checks for gaps in the frequency channels and that all MSs have the same frequency axis
        """
        # check that all MSs have the same frequency axis
        for MS_id in xrange(1,self.numMS):
            sw = pt.table(self.files[MS_id]+'::SPECTRAL_WINDOW', ack=False)
            if self.freq != sw.col('REF_FREQUENCY')[0] or self.nchan != sw.col('NUM_CHAN')[0] \
                    or not np.array_equal(self.chan_freqs_hz, sw.getcell('CHAN_FREQ',0)) \
                    or not np.array_equal(self.chan_width_hz, sw.getcell('CHAN_WIDTH',0)[0] ):
                self.log.critical('Frequency axis for MS {0} differs from the one for MS {1}! '
                                  'Exiting!'.format(self.files[MS_id],self.files[0]))
                sys.exit(1)
            sw.close()

        # check for gaps in the frequency channels
        self.missing_channels = []
        for i, (freq1, freq2) in enumerate(zip(self.chan_freqs_hz[:-1], self.chan_freqs_hz[1:])):
            ngap = int(round((freq2 - freq1)/self.chan_width_hz))
            self.missing_channels.extend([i + j + 1 for j in range(ngap-1)])
        self.log.debug('Missing channels: {}'.format(self.missing_channels))


    def chunk_input_files(self, chunksize, dirindparmdb, local_dir=None,
        clobber=True, test_run=False):
        """
        Make copies of input files that are smaller than 2*chunksize

        Chops off chunk of chunksize length until remainder is smaller than 2*chunksize
        Generates new self.files, self.msnames, and self.dirindparmdbs
        The direction independent parmDBs are fully copied into the new MSs

        Parameters
        ----------
        chunksize : float
            length of a chunk in seconds
        dirindparmdb : str
            Name of direction-independent instrument parmdb inside the new chunk files
        local_dir : str
            Path to local scratch directory for temp output. The file is then
            copied to the original output directory
        clobber : bool, optional
            If True, remove existing files if needed.
        test_run : bool, optional
            If True, don't actually do the chopping.

        """
        newfiles = []
        newdirindparmdbs = []
        for MS_id in xrange(self.numMS):
            nchunks = 1
            tab = pt.table(self.files[MS_id], ack=False)

            # Make filter for data columns that we don't need. These include imaging
            # columns and those made during initial subtraction
            colnames = tab.colnames()
            colnames_to_remove = ['MODEL_DATA', 'CORRECTED_DATA', 'IMAGING_WEIGHT',
                'SUBTRACTED_DATA_HIGH', 'SUBTRACTED_DATA_ALL_NEW', 'SUBTRACTED_DATA']
            colnames_to_keep = [c for c in colnames if c not in colnames_to_remove]

            timepersample = tab.getcell('EXPOSURE',0)
            timetab = tab.sort('unique desc TIME')
            tab.close()
            timearray = timetab.getcol('TIME')
            timetab.close()
            numsamples = len(timearray)
            mystarttime = np.min(timearray)
            myendtime = np.max(timearray)
            assert (timepersample*(numsamples-1)+.5) > (myendtime-mystarttime)
            if (myendtime-mystarttime) > (2.*chunksize):
                nchunks = int((numsamples*timepersample)/chunksize)
            if test_run:
                self.log.debug('Would split (or not) {0} into {1} chunks. '.format(self.files[MS_id], nchunks))
                tab.close()
                continue

            # Define directory where chunks are stored
            newdirname = self.chunks_dir
            if not os.path.exists(newdirname):
                os.mkdir(newdirname)

            if nchunks > 1:
                self.log.debug('Spliting {0} into {1} chunks...'.format(self.files[MS_id], nchunks))

                pool = multiprocessing.Pool()
                results = pool.map(process_chunk_star,
                    itertools.izip(itertools.repeat(self.files[MS_id]),
                    itertools.repeat(self.dirindparmdbs[MS_id]),
                    range(nchunks), itertools.repeat(nchunks),
                    itertools.repeat(mystarttime),
                    itertools.repeat(myendtime), itertools.repeat(chunksize),
                    itertools.repeat(dirindparmdb),
                    itertools.repeat(colnames_to_keep),
                    itertools.repeat(newdirname),
                    itertools.repeat(local_dir), itertools.repeat(clobber)))
                pool.close()
                pool.join()

                for chunk_file, chunk_parmdb in results:
                    newfiles.append(chunk_file)
                    newdirindparmdbs.append(chunk_parmdb)
            else:
                # Just copy the original files
                chunk_name = '{0}_chunk0.ms'.format(os.path.splitext(os.path.basename(self.files[MS_id]))[0])
                chunk_file = os.path.join(newdirname, chunk_name)
                if not os.path.exists(chunk_file):
                    shutil.copytree(self.files[MS_id], chunk_file)

                newdirindparmdb = os.path.join(chunk_file, dirindparmdb)
                if not os.path.exists(newdirindparmdb):
                    shutil.copytree(self.dirindparmdbs[MS_id], newdirindparmdb)

                newfiles.append(chunk_file)
                newdirindparmdbs.append(newdirindparmdb)

        # Check that each file has at least 10% unflagged data. If not, remove
        # it from the file list
        min_fraction = 0.1
        for f, p in zip(newfiles[:], newdirindparmdbs[:]):
            if self.find_unflagged_fraction(f) < min_fraction:
                newfiles.remove(f)
                newdirindparmdbs.remove(p)
                self.log.debug('Skipping file {0} in further processing '
                    '(unflagged fraction < {1}%)'.format(f, min_fraction*100.0))

        if test_run:
            return
        self.files = newfiles
        self.msnames = [ os.path.basename(MS) for MS in self.files ]
        self.dirindparmdbs = newdirindparmdbs
        self.numMS = len(self.files)


    def get_nearest_frequstep(self, freqstep):
        """
        Gets the nearest frequstep

        Parameters
        ----------
        freqstep : int
            Target frequency step

        Returns
        -------
        optimum_step : int
            Optimum frequency step nearest to target step

        """
        # first generate a list of possible values for freqstep
        if not hasattr(self, 'freq_divisors'):
            tmp_divisors = []
            for step in range(self.nchan,0,-1):
                if (self.nchan % step) == 0:
                    tmp_divisors.append(step)
            self.freq_divisors = np.array(tmp_divisors)
        idx = np.argmin(np.abs(self.freq_divisors-freqstep))
        return self.freq_divisors[idx]


    def find_unflagged_fraction(self, ms_file):
        """
        Finds the fraction of data that is unflagged

        Parameters
        ----------
        ms_file : str
            Filename of input MS

        Returns
        -------
        unflagged_fraction : float
            Fraction of unflagged data

        """
        unflagged_fraction = pt.taql('CALC sum([select nfalse(FLAG) from {0}]) '
            '/ sum([select nelements(FLAG) from {0}])'.format(ms_file))[0]

        return unflagged_fraction



def process_chunk_star(inputs):
    """
    Simple helper function for pool.map
    """
    return process_chunk(*inputs)


def process_chunk(ms_file, ms_parmdb, chunkid, nchunks, mystarttime, myendtime, chunksize, dirindparmdb,
    colnames_to_keep, newdirname, local_dir=None, clobber=True):
    """
    Processes one time chunk of input ms_file and returns new file names

    Parameters
    ----------
    ms_file : str
        Input MS file to chunk
    ms_parmdb : str
        Input dir-independent parmdb for input MS file
    chunkid : int
        ID of chunk
    nchunks : int
        Total number of chunks
    mystarttime : float
        Start time of MS file
    myendtime : float
        End time of MS file
    chunksize : float
        length of a chunk in seconds
    dirindparmdb : str
        Name of direction-independent instrument parmdb inside the new chunk files
    colnames_to_keep : list
        List of column names to keep in output chunk
    newdirname : str
        Name of output directory
    local_dir : str
        Path to local scratch directory for temp output. The file is then
        copied to the original output directory
    clobber : bool, optional
        If True, remove existing files if needed.

    Returns
    -------
    chunk_file : str
        Filename of chunk MS.
    newdirindparmdb : str
        Filename of direction-independent instrument parmdb for chunk_file

    """
    chunk_name = '{0}_chunk{1}.ms'.format(os.path.splitext(os.path.basename(ms_file))[0], chunkid)
    chunk_file = os.path.join(newdirname, chunk_name)
    old_chunk_file = os.path.join(os.path.dirname(ms_file), 'chunks', chunk_name)

    starttime = mystarttime+chunkid*chunksize
    endtime = mystarttime+(chunkid+1)*chunksize
    if chunkid == 0:
        starttime -= chunksize
    if chunkid == (nchunks-1):
        endtime += 2.*chunksize
    tab = pt.table(ms_file, ack=False)
    seltab = tab.query('TIME >= ' + str(starttime) + ' && TIME < ' + str(endtime),
        sortlist='TIME,ANTENNA1,ANTENNA2', columns=','.join(colnames_to_keep))

    if os.path.exists(chunk_file):
        try:
            newtab = pt.table(chunk_file, ack=False)
            if len(newtab) == len(seltab):
                copy = False
            newtab.close()
        except:
            copy = True
            os.shutil.rmtree(chunk_file)
    elif os.path.exists(old_chunk_file):
        # For compatibility, also search in old location
        try:
            newtab = pt.table(old_chunk_file, ack=False)
            if len(newtab) == len(seltab):
                copy = False
                chunk_file = old_chunk_file
            newtab.close()
        except:
            copy = True
    else:
        copy = True

    newdirindparmdb = os.path.join(chunk_file, dirindparmdb)

    if copy:
        if local_dir is not None:
            # Set output to temp directory
            chunk_file_original = chunk_file
            chunk_file = os.path.join(local_dir, os.path.basename(chunk_file_original))
            if os.path.exists(chunk_file):
                shutil.rmtree(chunk_file)

        seltab.copy(chunk_file, True)

        if local_dir is not None:
            # Copy temp file to original output location and clean up
            chunk_file_destination_dir = os.path.dirname(chunk_file_original)
            os.system('/usr/bin/rsync -a {0} {1}'.format(chunk_file, chunk_file_destination_dir))
            if not os.path.samefile(chunk_file, chunk_file_original):
                shutil.rmtree(chunk_file)
            chunk_file = chunk_file_original

        shutil.copytree(ms_parmdb, newdirindparmdb)
    seltab.close()
    tab.close()

    return (chunk_file, newdirindparmdb)
