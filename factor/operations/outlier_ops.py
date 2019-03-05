"""
Module that holds all operations for outlier sources

Classes
-------
OutlierPeel : Operation
    Runs the calibration for peeling an outlier source. Must be run in series
    as writes are made to original datasets

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
    Operation to peel an outlier direction
    """
    def __init__(self, parset, bands, direction, name='OutlierPeel'):
        super(OutlierPeel, self).__init__(parset, bands, direction,
            name=name)

        # Define extra parameters needed for this operation
        self.direction.set_imcal_parameters(parset, bands)
        ms_files = [band.files for band in self.bands]
        ms_files_single = []
        for bandfiles in ms_files:
            for filename in bandfiles:
                ms_files_single.append(filename)
        skymodels = [band.skymodel_dirindep for band in self.bands]

        # Parset for slow gain solve
        if self.direction.solve_all_correlations:
            selfcal_caltype = 'fulljones'
            fourpol = True # plot all correlations
        else:
            selfcal_caltype = 'diagonal'
            fourpol = False

        # Task for smoothing
        if self.parset['calibration_specific']['spline_smooth2d']:
            smooth_amps_task = 'smooth_amps_spline'
        else:
            smooth_amps_task = 'smooth_amps'

        # Set mapfile for selfcal operations
        if self.local_selfcal_scratch_dir is not None:
            self.direction.concat_data_mapfile = 'make_concat_data_sync_mapfile.output.mapfile'
        else:
            self.direction.concat_data_mapfile = 'concat_data.output.mapfile'

        self.parms_dict.update({'ms_files_single': ms_files_single,
                                'ms_files_grouped' : str(ms_files),
                                'skymodels': skymodels,
                                'fourpol': fourpol,
                                'selfcal_caltype': selfcal_caltype,
                                'smooth_amps_task': smooth_amps_task})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object for later reference
        self.direction.input_files_single_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'input_files_single.mapfile')
        self.direction.verify_subtract_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'verify_subtract.break.mapfile')
        self.direction.dir_dep_parmdb_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'merge_normalized_selfcal_parmdbs.mapfile')
        self.direction.converted_parmdb_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'convert_normalized_merged_selfcal_parmdbs.mapfile')
        self.direction.dir_indep_skymodels_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'full_skymodels.mapfile')
        self.direction.selfcal_plots_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'make_selfcal_plots.mapfile')
        if self.direction.create_preapply_h5parm:
            self.direction.preapply_parmdb_mapfile = os.path.join(self.pipeline_mapfile_dir,
                'create_preapply_h5parm.mapfile')
        self.direction.sourcedb_new_facet_sources = os.path.join(self.pipeline_mapfile_dir,
            'make_sourcedb_new_facet_sources_for_facet_imaging.mapfile')
        self.direction.diff_models_field_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'predict_and_difference_models.mapfile')

        # Store results of verify_subtract check. This will work if the verification
        # was done using multiple bands although we use only one at the moment
        if (os.path.exists(self.direction.verify_subtract_mapfile) and not
            self.parset['calibration_specific']['skip_selfcal_check']):
            ok_mapfile = DataMap.load(self.direction.verify_subtract_mapfile)
            ok_flags = [ast.literal_eval(item.file) for item in ok_mapfile]
            if all(ok_flags):
                self.direction.selfcal_ok = True
            else:
                self.direction.selfcal_ok = False
        elif self.parset['calibration_specific']['skip_selfcal_check']:
            self.direction.selfcal_ok = True
        else:
            self.direction.selfcal_ok = False

        # Delete temp data
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'shift_cal.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'apply_dir_dep.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'average_pre.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'average_post.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'sorted_groups.mapfile_groups')]
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()
        self.cleanup()
