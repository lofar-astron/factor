"""
Module that holds all facet operations

The facet calibration steps are divided into operations on the basis of whether
or not they can be run in parallel or in series.

Classes
-------
FacetSelfcal : Operation
    Runs the selfcal and imaging of a facet. May be run in parallel
FacetSub : Operation
    Subtracts all facet sources from data. Must be run in series as writes are
    made to original datasets
FacetImage : Operation
    Images the entire facet. May be run in parallel

"""
import os
import ast
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap


class FacetSelfcal(Operation):
    """
    Operation to selfcal a direction

    Two pipelines can be run, depending on whether casapy or wsclean is used
    for the full facet imaging:

    facetselfcal_casa_pipeline.parset - runs selfcal with casa

    facetselfcal_pipeline.parset - runs selfcal with wsclean

    """
    def __init__(self, parset, bands, direction):
        super(FacetSelfcal, self).__init__(parset, bands, direction,
            name='FacetSelfcal')

        # Set the pipeline parset to use
        if self.parset['facet_imager'].lower() == 'casa' or self.parset['facet_imager'].lower() == 'casapy':
            # Set parset template to CASA parset
            self.pipeline_parset_template = '{0}_casa_pipeline.parset'.format(self.name)
        else:
            # Set parset template to default (i.e., WSClean) parset
            self.pipeline_parset_template = '{0}_pipeline.parset'.format(self.name)

        # Define extra parameters needed for this operation (beyond those
        # defined in the master Operation class and as attributes of the
        # direction object)
        ms_files = [band.file for band in self.bands]
        skymodels = [band.skymodel_dirindep for band in self.bands]
        dir_indep_parmdbs = [band.dirindparmdb for band in self.bands]
        loopcount = max(1, self.parset['max_selfcal_loops'])
        self.parms_dict.update({'ms_files': ms_files,
                                'skymodels': skymodels,
                                'dir_indep_parmdbs': dir_indep_parmdbs,
                                'loopcount': loopcount})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object for later use
        self.direction.input_bands_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'input_bands.mapfile')
        self.direction.shifted_model_data_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'corrupt_final_model.mapfile')
        self.direction.diff_models_field_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'shift_diff_model_to_field.mapfile')
        self.direction.dir_indep_parmdbs_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'dir_indep_instrument_parmdbs.mapfile')
        self.direction.dir_indep_skymodels_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'full_skymodels.mapfile')
        self.direction.dir_indep_facet_skymodels_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'make_facet_skymodels_all.mapfile')
        self.direction.dir_dep_parmdb_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'merge_selfcal_parmdbs.mapfile')
        self.direction.facet_image_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'final_image.mapfile')
        self.direction.facet_model_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'final_model_rootnames.mapfile')
        self.wsclean_modelimg_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'expand_wsclean_model')
        self.wsclean_modelimg_size_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'expand_wsclean_model_size')
        self.verify_subtract_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'verify_subtract.mapfile')

        # Store results of verify_subtract check. This will work if the verification
        # was done using multiple bands although we use only one at the moment
        if os.path.exists(self.verify_subtract_mapfile) and not self.parset['skip_selfcal_check']:
            ok_mapfile = DataMap.load(self.verify_subtract_mapfile)
            ok_flags = [ast.literal_eval(item.file) for item in ok_mapfile]
            if all(ok_flags):
                self.direction.selfcal_ok = True
            else:
                self.direction.selfcal_ok = False
        elif self.parset['skip_selfcal_check']:
            self.direction.selfcal_ok = True
        else:
            self.direction.selfcal_ok = False

        # Delete all data used only for selfcal as they're no longer needed.
        # Note: we keep the data if selfcal failed verification, so that the user
        # can check them for problems
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'corrupt_all_model_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'corrupt_final_model.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'predict_all_facet_sources.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'shift_cal.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_corr.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_blavg_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_averaged.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'chunk_files.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat0_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat1_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat2_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat3_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat4_input.mapfile')]
        if not self.parset['keep_avg_facet_data'] and self.direction.name != 'target':
            # Add averaged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to reimage by hand (e.g.,
            # with a different weighting). They are always kept for the target
            self.direction.cleanup_mapfiles.append(
                os.path.join(self.pipeline_mapfile_dir, 'concat_averaged_input.mapfile'))
        if not self.parset['keep_unavg_facet_data']:
            # Add unaveraged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to phase shift them to
            # another direction (e.g., to combine several facets together before
            # imaging them all at once)
            self.direction.cleanup_mapfiles.append(
                os.path.join(self.pipeline_mapfile_dir, 'shift_empty.mapfile'))
        if self.direction.selfcal_ok or not self.parset['exit_on_selfcal_failure']:
            self.log.info('Cleaning up files (direction: {})'.format(self.direction.name))
            self.direction.cleanup()


class FacetSub(Operation):
    """
    Operation to subtract improved model
    """
    def __init__(self, parset, bands, direction):
        super(FacetSub, self).__init__(parset, bands, direction,
            name='FacetSub')


class FacetSubReset(Operation):
    """
    Operation to reset the subtraction of improved model
    """
    def __init__(self, parset, bands, direction):
        super(FacetSubReset, self).__init__(parset, bands, direction,
            name='FacetSubReset')


class FacetImage(Operation):
    """
    Operation to make the full facet image

    Four pipelines can be run, depending on whether casapy or wsclean is used
    for imaging and on whether an improved model (from selfcal) exists:

    facetimage_skymodel_casa_pipeline.parset - runs imaging with original
        skymodel using casa

    facetimage_imgmodel_casa_pipeline.parset - runs imaging with improved
        model using casa

    facetimage_skymodel_pipeline.parset - runs imaging with original
        skymodel using wsclean

    facetimage_imgmodel_pipeline.parset - runs imaging with improved
        model using wsclean

    """
    def __init__(self, parset, bands, direction):
        super(FacetImage, self).__init__(parset, bands, direction,
            name='FacetImage')

        # Set the pipeline parset to use
        if self.parset['facet_imager'].lower() == 'casa':
            # Set parset template to casa parset
            if not self.direction.selfcal_ok:
                # Set parset template to sky-model parset
                self.pipeline_parset_template = '{0}_skymodel_casa_pipeline.parset'.format(self.name)
            else:
                # Set parset template to facet model-image parset
                self.pipeline_parset_template = '{0}_imgmodel_casa_pipeline.parset'.format(self.name)
        else:
            # Set parset template to wsclean parset
            if not self.direction.selfcal_ok:
                # Set parset template to sky-model parset
                self.pipeline_parset_template = '{0}_skymodel_pipeline.parset'.format(self.name)
            else:
                # Set parset template to facet model-image parset
                self.pipeline_parset_template = '{0}_imgmodel_pipeline.parset'.format(self.name)

        # Define extra parameters needed for this operation (beyond those
        # defined in the master Operation class and as attributes of the
        # direction object)
        ms_files = [band.file for band in self.bands]
        skymodels = [band.skymodel_dirindep for band in self.bands]
        dir_indep_parmdbs = [band.dirindparmdb for band in self.bands]
        self.parms_dict.update({'ms_files': ms_files,
                                'dir_indep_parmdbs': dir_indep_parmdbs,
                                'skymodels': skymodels})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object for later use
        self.direction.facet_image_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'final_image.mapfile')

        # Delete temp data
        if not self.parset['keep_avg_facet_data'] and self.direction.name != 'target':
            # Add averaged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to reimage by hand (e.g.,
            # with a different weighting). They are always kept for the target
            self.direction.cleanup_mapfiles.append(
                os.path.join(self.pipeline_mapfile_dir, 'concat_averaged_input.mapfile'))
        if not self.parset['keep_unavg_facet_data']:
            # Add unaveraged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to phase shift them to
            # another direction (e.g., to combine several facets together before
            # imaging them all at once)
            self.direction.cleanup_mapfiles.append(
                os.path.join(self.pipeline_mapfile_dir, 'shift_empty.mapfile'))
        self.direction.cleanup()

