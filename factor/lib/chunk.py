"""
Definition of the time chunk class
"""
import logging
import os

log = logging.getLogger('parset')

class Chunk(object):
    """
    The Chunk object contains parameters needed for each time chunk (MS)
    """
    def __init__(self, MSfile, index):
        """
        Create Chunk object

        Parameters
        ----------
        MSfile : str
            Filename of MS
        index : int
            Index of chunk

        """
        self.parent_file = MSfile
        self.index = index
        self.file = '{0}-chunk_{1}.ms'.format(os.path.splitext(MSfile)[0], index)
        self.msname = self.file.split('/')[-1]
        self.name = 'chunk_{0}'.format(index)
        self.parmdb = 'instrument'
        self.parmdb_phaseamp1a = 'instrument_phase1'
        self.parmdb_phaseamp2a = 'instrument_phase2'
        self.parmdb_phaseamp1b = 'instrument_amp1'
        self.parmdb_phaseamp2b = 'instrument_amp2'
