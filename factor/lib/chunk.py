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
    def __init__(self, op_parset, MSfile, index, prefix=None, direction=None,
        outdir=None):
        """
        Create Chunk object

        Parameters
        ----------
        op_parset : str
            Parset of calling operation
        MSfile : str
            Filename of MS
        index : int
            Index of chunk
        prefix : str, optional
            Prefix to use for names
        direction : Direction object, optional
            Direction for this chunk
        outdir : str, optional
            Absolute path to output directory. If None, the directory of the
            parent is used

        """
        self.parent_file = MSfile
        self.index = index
        self.direction = direction
        if outdir is None:
            self.file = '{0}-chunk_{1}.ms'.format(os.path.splitext(MSfile)[0],
                index)
        else:
            self.file = os.path.join(outdir, '{0}-chunk_{1}.ms'.format(
                os.path.splitext(os.path.basename(MSfile))[0], index))
        self.msname = self.file.split('/')[-1]
        self.name = 'chunk_{0}'.format(index)

        factor_working_dir = op_parset['dir_working']
        op_name = op_parset['op_name']
        self.parmdb_dir = '{0}/parmdbs/{1}/'.format(factor_working_dir,
            op_name)
        if direction is not None:
            self.parmdb_dir += '{0}/'.format(direction.name)
        self.parmdb_dir += 'chunks/'
        if not os.path.exists(self.parmdb_dir):
            os.makedirs(self.parmdb_dir)

        self.parmdb_phaseonly1 = self.parmdb_dir + 'chunk{0}_instrument_phase1'.format(self.index)
        self.parmdb_phaseonly2 = self.parmdb_dir + 'chunk{0}_instrument_phase2'.format(self.index)
        self.parmdb_phaseamp_phase1 = self.parmdb_dir + 'chunk{0}_instrument_phaseamp_phase1'.format(self.index)
        self.parmdb_phaseamp_phase2 = self.parmdb_dir + 'chunk{0}_instrument_phaseamp_phase2'.format(self.index)
        self.parmdb_phaseamp_amp1 = self.parmdb_dir + 'chunk{0}_instrument_phaseamp_amp1'.format(self.index)
        self.parmdb_phaseamp_amp2 = self.parmdb_dir + 'chunk{0}_instrument_phaseamp_amp2'.format(self.index)

    def copy_from_parent(column_list):
        """
        Copy columns from the parent MS file to the chunk file

        The columns in both parent and chunk file will have the same names
        """
        from factor.lib.operation_lib import copy_column

        for column in column_list:
            copy_column(self.parent_file, column, column, ms_from=self.file)


    def copy_to_parent(column_list):
        """
        Copy columns from the chunk file to the parent MS file

        The columns in both parent and chunk file will have the same names
        """
        from factor.lib.operation_lib import copy_column

        for column in column_list:
            copy_column(self.file, column, column, ms_from=self.parent_file)


