"""
Definition of the band class
"""
import os
import sys
import shutil
import logging
import pyrap.tables as pt
import lofar.parmdb
import numpy as np


class Band(object):
    """
    The Band object contains parameters needed for each band (MS)

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
    test_run : bool, optional
        If True, use test image sizes

    """
    def __init__(self, MSfiles, factor_working_dir, dirindparmdb,
        skymodel_dirindep=None, test_run=False):
        """
        Create Band object

        Parameters
        ----------
        MSfiles : list of str
            Filenames of MSs for the same frequency band
        factor_working_dir : str
            Full path of working directory
        dirindparmdb : str
            Name of direction-independent instrument parmdb (relative to MSfile)
        skymodel_dirindep : str
            Full path of direction-independent sky model
        test_run : bool, optional
            If True, use test image sizes

        """
        self.files = MSfiles
        self.msnames = [ MS.split('/')[-1] for MS in self.files ]
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
        self.log.debug('MS filename is {}'.format(self.msname))

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

        # cut input files into chunks if needed
        chunksize = 2400. # in seconds -> 40min
        self.chunk_input_files(chunksize, dirindparmdb, clobber=True, test_run=test_run)

        # Check for SUBTRACTED_DATA_ALL column and calculate times and number
        # of samples
        self.has_sub_data = True
        self.sumsamples = 0
        self.minSamplesPerFile = 4294967295  # If LOFAR lasts that many seconds then I buy you a beer.
        self.starttime = np.finfo('d').max
        self.endtime = 0.
        for MSid in xrange(self.numMS):
            tab = pt.table(self.files[MSid], ack=False)
            if not 'SUBTRACTED_DATA_ALL' in tab.colnames():
                self.has_sub_data = False
            self.has_sub_data_new = False
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

        # Set the initsubtract cell sizes to default / dummy values
        self.cellsize_highres_deg = 0.00208 # initsubtract high-res cell size
        self.cellsize_lowres_deg = 0.00694 # initsubtract low-res cell size

        self.completed_operations = []
        self.skip = False
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


    def convert_parmdb_to_phasors_id(self,pdb_id=0):
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
        solnames = pdb.getNames()
        for name in solnames:
            if solname[0:9] != 'Gain:0:0:' and solname[0:9] != 'Gain:1:1:':
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



    def chunk_input_files(self, chunksize, dirindparmdb, clobber=True, test_run=False):
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
        clobber : bool, optional
            If True, remove existing files if needed.      
        test_run : bool, optional
            If True, don't actually do the choping.
        """
        newfiles = []
        newdirindparmdbs = []
        for MS_id in xrange(self.numMS):
            nchunks = 1
            tab = pt.table(self.files[MS_id], ack=False)            
            timepersample = tab.getcell('EXPOSURE',0)
            timetab = tab.sort('unique desc TIME')
            timearray = timetab.getcol('TIME')
            numsamples = len(timearray)
            mystarttime = np.min(timearray)
            myendtime = np.max(timearray)
            assert (timepersample*(numsamples-1)+.5) > (myendtime-mystarttime)
            if (myendtime-mystarttime) > (2.*chunksize):
                nchunks = int((numsamples*timepersample)/chunksize)
            if test_run:
                self.log.debug('Would split (or not) {0} into {1} chunks. '.format(self.files[MS_id],nchunks))
                tab.close()
                continue
            if nchunks > 1:
                newdirname = os.path.join(os.path.dirname(self.files[MS_id]),'chunks')
                if not os.path.exists(newdirname):
                    os.mkdir(newdirname)
                for chunkid in range(nchunks):
                    chunk_name = '{0}_chunk{1}.ms'.format(os.path.splitext(os.path.basename(self.files[MS_id]))[0], chunkid)
                    chunk_file = os.path.join(newdirname,chunk_name)
                    newdirindparmdb = os.path.join(chunk_file, dirindparmdb)
                    starttime = mystarttime+chunkid*chunksize
                    endtime = mystarttime+(chunkid+1)*chunksize
                    if chunkid == 0: 
                        starttime -= chunksize
                    if chunkid == (nchunks-1):
                        endtime += 2.*chunksize
                    seltab = tab.query('TIME >= ' + str(starttime) + ' && '
                                       'TIME < ' + str(endtime), sortlist='TIME,ANTENNA1,ANTENNA2')
                    self.log.debug('Going to copy {0} samples to file {1}'.format(str(len(seltab)),chunk_file))
                    if os.path.exists(chunk_file):
                        try:
                            newtab = pt.table(chunk_file, ack=False)
                            if len(newtab) == len(seltab):
                                self.log.debug('Found existing file of correct length, not copying!')
                                copy = False
                            newtab.close()
                        except:
                            copy = True
                            os.shutil.rmtree(chunk_file)
                    else:
                        copy = True
                    if copy:
                        seltab.copy(chunk_file, True)
                        shutil.copytree(self.dirindparmdbs[MS_id],newdirindparmdb)
                    seltab.close()
                    newfiles.append(chunk_file)
                    newdirindparmdbs.append(newdirindparmdb)
            else:
               newfiles.append(self.files[MS_id])
               newdirindparmdbs.append(self.dirindparmdbs[MS_id])
            tab.close()            
        if test_run:
            return
        self.files = newfiles
        self.msnames = [ MS.split('/')[-1] for MS in self.files ]
        self.dirindparmdbs = newdirindparmdbs
        self.numMS = len(self.files)

  

    def set_image_sizes(self, test_run=False,cellsize_highres_deg=None,cellsize_lowres_deg=None,
                        fieldsize_highres=2.5,fieldsize_lowres=6.5):
        """
        Sets sizes for initsubtract images

        The image sizes are scaled from the mean primary-beam FWHM. For
        the high-res image, we use 2.5 * FWHM; for low-res, we use 6.5 * FHWM.

        Parameters
        ----------
        test_run : bool, optional
            If True, use test sizes
        cellsize_highres_deg : float, optional
            cellsize for the high-res images in deg
        cellsize_lowres_deg : float, optional
            cellsize for the low-res images in deg
        fieldsize_highres : float, optional
            How many FWHM's shall the high-res images be.
        fieldsize_lowres : float, optional
            How many FWHM's shall the low-res images be.
        """
        if cellsize_highres_deg:
            self.cellsize_highres_deg = cellsize_highres_deg
        if cellsize_lowres_deg:
            self.cellsize_lowres_deg = cellsize_lowres_deg
        if not test_run:
            if not hasattr(self, 'mean_el_rad'):
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

            # Calculate mean FOV
            sec_el = 1.0 / np.sin(self.mean_el_rad)
            self.fwhm_deg = 1.1 * ((3.0e8 / self.freq) / self.diam) * 180. / np.pi * sec_el
            self.imsize_high_res = self.get_optimum_size(self.fwhm_deg
                /self.cellsize_highres_deg * fieldsize_highres)
            self.imsize_low_res = self.get_optimum_size(self.fwhm_deg
                /self.cellsize_lowres_deg * fieldsize_lowres)
        else:
            self.imsize_high_res = self.get_optimum_size(128)
            self.imsize_low_res = self.get_optimum_size(128)


    def get_optimum_size(self, size):
        """
        Gets the nearest optimum image size

        Taken from the casa source code (cleanhelper.py)

        Parameters
        ----------
        size : int
            Target image size in pixels

        Returns
        -------
        optimum_size : int
            Optimum image size nearest to target size

        """
        import numpy

        def prime_factors(n, douniq=True):
            """ Return the prime factors of the given number. """
            factors = []
            lastresult = n
            sqlast=int(numpy.sqrt(n))+1
            if n == 1:
                return [1]
            c=2
            while 1:
                 if (lastresult == 1) or (c > sqlast):
                     break
                 sqlast=int(numpy.sqrt(lastresult))+1
                 while 1:
                     if(c > sqlast):
                         c=lastresult
                         break
                     if lastresult % c == 0:
                         break
                     c += 1

                 factors.append(c)
                 lastresult /= c

            if (factors==[]): factors=[n]
            return  numpy.unique(factors).tolist() if douniq else factors

        n = int(size)
        if (n%2 != 0):
            n+=1
        fac=prime_factors(n, False)
        for k in range(len(fac)):
            if (fac[k] > 7):
                val=fac[k]
                while (numpy.max(prime_factors(val)) > 7):
                    val +=1
                fac[k]=val
        newlarge=numpy.product(fac)
        for k in range(n, newlarge, 2):
            if ((numpy.max(prime_factors(k)) < 8)):
                return k
        return newlarge


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
