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
        self.parmdb = '{0}/instrument'.format(self.file)
        self.parmdb_phaseamp1a = self.parmdb + '_phase1'
        self.parmdb_phaseamp2a = self.parmdb + '_phase2'
        self.parmdb_phaseamp1b = self.parmdb + '_amp1'
        self.parmdb_phaseamp2b = self.parmdb + '_amp2'
