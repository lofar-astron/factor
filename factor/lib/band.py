"""
Definition of the band class
"""
import os
import sys
import logging
import pyrap.tables as pt
import lofar.parmdb
import numpy as np


class Band(object):
    """
    The Band object contains parameters needed for each band (MS)
    """
    def __init__(self, MSfile, factor_working_dir, dirindparmdb,
        skymodel_dirindep=None, test_run=False):
        """
        Create Band object

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
        self.file = MSfile
        self.msname = self.file.split('/')[-1]
        self.dirindparmdb = os.path.join(self.file, dirindparmdb)
        self.skymodel_dirindep = skymodel_dirindep

        # Get the frequency info and set name
        sw = pt.table(self.file+'::SPECTRAL_WINDOW', ack=False)
        self.freq = sw.col('REF_FREQUENCY')[0]
        self.nchan = sw.col('NUM_CHAN')[0]
        self.chan_freqs_hz = sw.col('CHAN_FREQ')[0]
        self.chan_width_hz = sw.col('CHAN_WIDTH')[0][0]
        sw.close()
        self.name = 'Band_{0:.2f}MHz'.format(self.freq/1e6)
        self.log = logging.getLogger('factor.{}'.format(self.name))

        # Do some checks
        self.check_freqs()
        self.check_parmdb()

        # Get the field RA and Dec
        obs = pt.table(self.file+'::FIELD', ack=False)
        self.ra = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][0]))
        if self.ra < 0.:
            self.ra = 360.0 + (self.ra)
        self.dec = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][1]))
        obs.close()

        # Get the station diameter
        ant = pt.table(self.file+'::ANTENNA', ack=False)
        self.diam = float(ant.col('DISH_DIAMETER')[0])
        ant.close()

        # Check for SUBTRACTED_DATA_ALL column and calculate times and number
        # of samples
        tab = pt.table(self.file, ack=False)
        if 'SUBTRACTED_DATA_ALL' in tab.colnames():
            self.has_sub_data = True
        else:
            self.has_sub_data = False
        self.has_sub_data_new = False
        self.starttime = tab.col('TIME')[0]
        self.endtime = tab.col('TIME')[-1]
        for t2 in tab.iter(["ANTENNA1","ANTENNA2"]):
            if (t2.getcell('ANTENNA1',0)) < (t2.getcell('ANTENNA2',0)):
                self.timepersample = t2.col('TIME')[1] - t2.col('TIME')[0]
                self.nsamples = t2.nrows()
                break
        tab.close()

        # Set the initsubtract cell sizes
        self.cellsize_highres_deg = 0.00208 # initsubtract high-res cell size
        self.cellsize_lowres_deg = 0.00694 # initsubtract low-res cell size

        self.completed_operations = []
        self.skip = False


    def check_parmdb(self):
        """
        Checks the dir-indep instrument parmdb for various problems
        """
        # Check for special BBS table name "instrument"
        if os.path.basename(self.dirindparmdb) == 'instrument':
            self.dirindparmdb += '_dirindep'
            if not os.path.exists(self.dirindparmdb):
                if not os.path.exists(os.path.join(self.file, 'instrument')):
                    self.log.critical('Direction-independent instument parmdb not found '
                        'for band {0}'.format(self.file))
                    sys.exit(1)
                self.log.warn('Direction-independent instument parmdb for band {0} is '
                    'named "instrument". Copying to "instrument_dirindep" so that BBS '
                    'will not overwrite this table...'.format(self.file))
                os.system('cp -r {0} {1}'.format(os.path.join(self.file,
                    'instrument'), self.dirindparmdb))
        if not os.path.exists(self.dirindparmdb):
            self.log.critical('Direction-independent instrument parmdb "{0}" not found '
                'for band {1}'.format(self.dirindparmdb, self.file))
            sys.exit(1)

        # Check whether there are ampl/phase or real/imag
        try:
            pdb = lofar.parmdb.parmdb(self.dirindparmdb)
            solname = pdb.getNames()[0]
        except IndexError:
            self.log.critical('Direction-independent instument parmdb appears to be empty '
                        'for band {0}'.format(self.file))
            sys.exit(1)
        if 'Real' in solname or 'Imag' in solname:
            # Convert real/imag to phasors
            self.log.warn('Direction-independent instument parmdb for band {0} contains '
                'real/imaginary values. Converting to phase/amplitude...'.format(self.file))
            self.convert_parmdb_to_phasors()

        # Check that there aren't extra default values in the parmdb, as this
        # confuses DPPP
        defvals = pdb.getDefValues()
        for v in defvals:
            if 'Ampl' not in v and 'Phase' not in v:
                pdb.deleteDefValues(v)
        pdb.flush()


    def convert_parmdb_to_phasors(self):
        """
        Converts instrument parmdb from real/imag to phasors
        """

        phasors_parmdb_file = self.dirindparmdb + '_phasors'
        if os.path.exists(phasors_parmdb_file):
            return

        pdb_in = lofar.parmdb.parmdb(self.dirindparmdb)
        pdb_out = lofar.parmdb.parmdb(phasors_parmdb_file, create=True)

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
        self.dirindparmdb = phasors_parmdb_file


    def check_freqs(self):
        """
        Checks for gaps in the frequency channels
        """
        self.missing_channels = []
        for i, (freq1, freq2) in enumerate(zip(self.chan_freqs_hz[:-1], self.chan_freqs_hz[1:])):
            ngap = int(round((freq2 - freq1)/self.chan_width_hz))
            self.missing_channels.extend([i + j + 1 for j in range(ngap)])


    def set_image_sizes(self, test_run=False):
        """
        Sets sizes for initsubtract images

        The image sizes are scaled from the mean primary-beam FWHM. For
        the high-res image, we use 2.5 * FWHM; for low-res, we use 6.5 * FHWM.

        Parameters
        ----------
        test_run : bool, optional
            If True, use test sizes

        """
        if not test_run:
            if not hasattr(self, 'mean_el_rad'):
                # Add (virtual) elevation column to MS
                try:
                    pt.addDerivedMSCal(self.file)
                except RuntimeError:
                    # RuntimeError indicates column already exists
                    pass

                # Check for SUBTRACTED_DATA_ALL column and calculate mean elevation
                tab = pt.table(self.file, ack=False)
                self.mean_el_rad = np.mean(tab.getcol('AZEL1', rowincr=10000)[:, 1])
                tab.close()

                # Remove (virtual) elevation column from MS
                pt.removeDerivedMSCal(self.file)

            # Calculate mean FOV
            sec_el = 1.0 / np.sin(self.mean_el_rad)
            self.fwhm_deg = 1.1 * ((3.0e8 / self.freq) / self.diam) * 180. / np.pi * sec_el
            self.imsize_high_res = self.get_optimum_size(self.fwhm_deg
                /self.cellsize_highres_deg* 2.5)
            self.imsize_low_res = self.get_optimum_size(self.fwhm_deg
                /self.cellsize_lowres_deg*6.5)
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
