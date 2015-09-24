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
        if self.parset['facet_imager'].lower() == 'casa':
            # Set parset template to casa parset
            self.pipeline_parset_template = '{0}_casa_pipeline.parset'.format(self.name)
        else:
            # Set parset template to wsclean parset
            self.pipeline_parset_template = '{0}_pipeline.parset'.format(self.name)

        # Define extra parameters needed for this operation (beyond those
        # defined in the master Operation class and as attributes of the
        # direction object)
        ms_files = [band.file for band in self.bands]
        skymodels = [band.skymodel_dirindep for band in self.bands]
        dir_indep_parmdbs = [band.dirindparmdb for band in self.bands]
        if self.direction.use_new_sub_data:
            subtracted_data_colname = 'SUBTRACTED_DATA_ALL_NEW'
        else:
            subtracted_data_colname = 'SUBTRACTED_DATA_ALL'
        if self.direction.nchannels > 1:
            nterms = 2
            casa_suffix = '.tt0'
            wsclean_suffix = '-MFS-image.fits'
        else:
            nterms = 1
            casa_suffix = None
            wsclean_suffix = '-image.fits'
        self.parms_dict.update({'ms_files': ms_files,
                                'subtracted_data_colname': subtracted_data_colname,
                                'skymodels': skymodels,
                                'dir_indep_parmdbs': dir_indep_parmdbs,
                                'casa_suffix': casa_suffix,
                                'wsclean_suffix': wsclean_suffix,
                                'nterms': nterms})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object for later use
        self.direction.input_bands_datamap = os.path.join(self.mapfile_dir,
            'input_bands.datamap')
        self.direction.shifted_bands_datamap = os.path.join(self.mapfile_dir,
            'shifted_empty_bands.datamap')
        self.direction.diff_models_field_datamap = os.path.join(self.mapfile_dir,
            'diff_models_field.datamap')
        self.direction.dir_indep_parmdbs_datamap = os.path.join(self.mapfile_dir,
            'dir_indep_instrument_parmdbs.datamap')
        self.direction.dir_indep_skymodels_datamap = os.path.join(self.mapfile_dir,
            'full_skymodels.datamap')
        self.direction.dir_indep_facet_skymodels_datamap = os.path.join(self.mapfile_dir,
            'facet_skymodels.datamap')
        self.direction.dir_dep_parmdb_datamap = os.path.join(self.mapfile_dir,
            'dir_dep_parmdb.datamap')
        self.direction.facet_image_mapfile = os.path.join(self.mapfile_dir,
            'final_image.datamap')
        self.direction.facet_model_mapfile = os.path.join(self.mapfile_dir,
            'final_model_rootnames.datamap')
        self.verify_subtract_OK_mapfile = os.path.join(self.mapfile_dir,
            'verify_subtract_OK.datamap')

        # Store results of verify_subtract check
        if os.path.exists(self.verify_subtract_OK_mapfile):
            ok_datamap = DataMap.load(self.verify_subtract_OK_mapfile)
            ok_flags = [ast.literal_eval(item.file) for item in ok_datamap]
            if all(ok_flags):
                self.direction.selfcal_ok = True
            else:
                self.direction.selfcal_ok = False
        else:
            self.direction.selfcal_ok = False

        # Delete all data used only for selfcal as they're no longer needed.
        # Note: we keep the data if selfcal failed verification, so that the user
        # can check them for problems. The data used to image the full facet
        # are also kept in case the user wants to re-image the facet by hand
        self.direction.cleanup_mapfiles = [
            os.path.join(self.mapfile_dir, 'shifted_cal_bands.datamap'),
            os.path.join(self.mapfile_dir, 'chunk_files.datamap'),
            os.path.join(self.mapfile_dir, 'concat_input.datamap'),
            os.path.join(self.mapfile_dir, 'concat0_input.datamap'),
            os.path.join(self.mapfile_dir, 'concat1_input.datamap'),
            os.path.join(self.mapfile_dir, 'concat2_input.datamap'),
            os.path.join(self.mapfile_dir, 'concat3_input.datamap'),
            os.path.join(self.mapfile_dir, 'concat4_input.datamap')]
        if self.direction.selfcal_ok:
            self.direction.cleanup()


class FacetSub(Operation):
    """
    Operation to subtract improved model

    Two pipelines can be run, depending on whether the empty datasets must
    be updated before subtraction is done:

    facetsub_single_pipeline.parset - runs the subraction without updating the
        empty datasets

    facetsub_pipeline.parset - runs the subtraction after updating the emtpy
        datasets

    """
    def __init__(self, parset, direction):
        super(FacetSub, self).__init__(parset, None, direction,
            name='FacetSub')


    def finalize(self):
        """
        Finalize this operation
        """
        # Add model datamap to direction object for later use
        self.direction.facet_model_data_mapfile = os.path.join(self.mapfile_dir,
            'shifted_to_field_models.datamap')

        # Delete shifted data as they're no longer needed
        self.direction.cleanup_mapfiles = [
            self.direction.shifted_bands_datamap,
            os.path.join(self.mapfile_dir, 'shifted_to_field_subdata.datamap')]
        self.direction.cleanup()


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
        if self.direction.nchannels > 1:
            nterms = 2
            casa_suffix = '.tt0'
            wsclean_suffix = '-MFS-image.fits'
        else:
            nterms = 1
            casa_suffix = None
            wsclean_suffix = '-image.fits'
        self.parms_dict.update({'ms_files': ms_files,
                                'dir_indep_parmdbs': dir_indep_parmdbs,
                                'skymodels': skymodels,
                                'casa_suffix': casa_suffix,
                                'wsclean_suffix': wsclean_suffix,
                                'nterms': nterms})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object for later use
        self.direction.facet_image_mapfile = os.path.join(self.mapfile_dir,
            'final_image.datamap')
        self.direction.facet_model_mapfile = os.path.join(self.mapfile_dir,
            'final_model_rootnames.datamap')

        # Delete shifted data as they're no longer needed
        self.direction.cleanup_mapfiles = [
            os.path.join(self.mapfile_dir, 'shifted_empty_final_bands.datamap'),
            os.path.join(self.mapfile_dir, 'shifted_all_final_bands.datamap'),
            os.path.join(self.mapfile_dir, 'shifted_avg_final_bands.datamap')]
        self.direction.cleanup()
