"""
Definition of the band class
"""
import os
import pyrap.tables as pt
import numpy as np


class Band(object):
    """
    The Band object contains parameters needed for each band (MS)
    """
    def __init__(self, MSfile, factor_working_dir, test_run=False):
        """
        Create Band object

        Parameters
        ----------
        MSfile : str
            Filename of MS
        factor_working_dir : str
            Full path of working directory
        test_run : bool, optional
            If True, use test image sizes

        """
        self.file = MSfile
        self.msname = self.file.split('/')[-1]

        # Get the reference frequency and set name
        sw = pt.table(self.file+'::SPECTRAL_WINDOW', ack=False)
        self.freq = sw.col('REF_FREQUENCY')[0]
        sw.close()
        self.name = 'Band_{0:.2f}MHz'.format(self.freq/1e6)

        # Get the field RA and Dec
        obs = pt.table(self.file+'::FIELD', ack=False)
        self.ra = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][0]))
        if self.ra < 0.:
            self.ra = 360.0 + (self.ra)
        self.dec = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][1]))
        obs.close()

        # Get the station diameter
        ant = pt.table(self.file+'::ANTENNA', ack=False)
        diam = float(ant.col('DISH_DIAMETER')[0])
        ant.close()

        # Add (virtual) elevation column to MS
        try:
            pt.addDerivedMSCal(self.file)
        except RuntimeError:
            pass

        # Check for SUBTRACTED_DATA_ALL column and calculate mean elevation
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
        self.mean_el_rad = np.mean(tab.getcol('AZEL1', rowincr=10000)[:, 1])
        tab.close()

        # Calculate mean FOV
        sec_el = 1.0 / np.sin(self.mean_el_rad)
        self.fwhm_deg = 1.1 * ((3.0e8 / self.freq) / diam) * 180. / np.pi * sec_el

        # Set the initsubtract image sizes from the primary-beam FWHM. For
        # high-res image, use 2.5 * FWHM; for low-res, use 6.5 * FHWM
        self.cellsize_highres_deg = 0.00208 # initsubtract high-res cell size
        self.cellsize_lowres_deg = 0.00694 # initsubtract low-res cell size
        if not test_run:
            self.imsize_high_res = self.get_optimum_size(self.fwhm_deg
                / self.cellsize_highres_deg*2.5)
            self.imsize_low_res = self.get_optimum_size(self.fwhm_deg
                / self.cellsize_lowres_deg*6.5)
        else:
            self.imsize_high_res = self.get_optimum_size(128)
            self.imsize_low_res = self.get_optimum_size(128)

        self.completed_operations = []
        self.save_file = os.path.join(factor_working_dir, 'state',
            self.name+'_save.pkl')


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


    def save_state(self):
        """
        Saves the state to a file
        """
        import pickle

        with open(self.save_file, 'wb') as f:
            pickle.dump(self.__dict__, f)


    def load_state(self):
        """
        Loads the state from a file

        Returns
        -------
        success : bool
            True if state was successfully loaded, False if not
        """
        import pickle

        try:
            with open(self.save_file, 'r') as f:
                self.__dict__ = pickle.load(f)
            return True
        except:
            return False
