"""
Module that holds all operations for outlier sources

The outlier calibration steps are divided into operations on the basis of whether
or not they can be run in parallel or in series.

Classes
-------
OutlierPeel : Operation
    Runs the calibration for peeling an outlier source. May be run in parallel
OutlierSub : Operation
    Subtracts outlier sources from data. Must be run in series as writes are
    made to original datasets

"""
import os
import ast
import sys
import logging
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap

log = logging.getLogger('factor:outlier_ops')

class OutlierPeel(Operation):
    """
    Operation to peel a direction
    """
    def __init__(self, parset, bands, direction):
        super(OutlierPeel, self).__init__(parset, bands, direction,
            name='OutlierPeel')

        # Define extra parameters needed for this operation (beyond those
        # defined in the master Operation class and as attributes of the
        # direction object)
        ms_files = [band.files for band in self.bands]
        ms_files_single = []
        for bandfiles in ms_files:
            for filename in bandfiles:
                ms_files_single.append(filename)
        dir_indep_parmDBs = []
        for band in self.bands:
            for parmdb in band.dirindparmdbs:
                dir_indep_parmDBs.append(parmdb)
        skymodels = [band.skymodel_dirindep for band in self.bands]

        self.parms_dict.update({'ms_files_single': ms_files_single,
                                'ms_files_grouped' : str(ms_files),
                                'skymodels': skymodels,
                                'dir_indep_parmDBs': dir_indep_parmDBs})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object for later reference
        self.direction.input_files_single_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'input_files_single.mapfile')
        self.direction.subtracted_data_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'add_all_facet_sources.mapfile')
        self.direction.verify_subtract_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'verify_subtract.break.mapfile')

        # Store results of verify_subtract check. This will work if the verification
        # was done using multiple bands although we use only one at the moment
        if os.path.exists(self.direction.verify_subtract_mapfile) and not self.parset['skip_selfcal_check']:
            ok_mapfile = DataMap.load(self.direction.verify_subtract_mapfile)
            ok_flags = [ast.literal_eval(item.file) for item in ok_mapfile]
            if all(ok_flags):
                self.direction.selfcal_ok = True
            else:
                self.direction.selfcal_ok = False
        elif self.parset['skip_selfcal_check']:
            self.direction.selfcal_ok = True
        else:
            self.direction.selfcal_ok = False

        # Delete temp data
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'shift_and_average.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_blavg_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'predict_outlier_model.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'corrupt_outlier_model.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'average_pre.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'average_post.mapfile')]
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()


class OutlierSub(Operation):
    """
    Operation to subtract improved model
    """
    def __init__(self, parset, bands, direction):
        super(OutlierSub, self).__init__(parset, bands, direction,
            name='OutlierSub')

        # Delete temp data
        self.direction.cleanup_mapfiles = [self.direction.subtracted_data_mapfile]
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()
