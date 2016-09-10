"""
Module that holds all operations for facets and patches

The facet calibration steps are divided into operations on the basis of whether
or not they can be run in parallel or in series.

Classes
-------
FacetSelfcal : Operation
    Runs the selfcal and imaging of a facet. May be run in parallel
FacetPeel : Operation
    Runs the peeling and imaging of a facet. May be run in parallel
FacetSub : Operation
    Subtracts all facet sources from data. Must be run in series as writes are
    made to original datasets
FacetSubreset : Operation
    Resets previous subtraction of all facet sources from data. Must be run in
    series as writes are made to original datasets
FacetImage : Operation
    Images the entire facet. May be run in parallel
FacetPeelImage:  Operation
    Images the entire facet after facet was peeled. May be run in parallel

"""
import os
import ast
from factor.lib.operation import Operation
from factor.operations.outlier_ops import OutlierPeel
from lofarpipe.support.data_map import DataMap


class FacetSelfcal(Operation):
    """
    Operation to selfcal a facet

    Two pipelines can be run, depending on whether multiresolution selfcal is
    desired or not:

    facetselfcal_pipeline.parset - normal selfcal

    facetselfcal_taper_pipeline.parset - multiresolution selfcal

    """
    def __init__(self, parset, bands, direction):
        super(FacetSelfcal, self).__init__(parset, bands, direction,
            name='FacetSelfcal')

        # Set the pipeline parset to use
        if self.parset['calibration_specific']['multires_selfcal']:
            # Set parset template to multi-resolution selfcal parset
            self.pipeline_parset_template = '{0}_taper_pipeline.parset'.format(self.name)
        else:
            self.pipeline_parset_template = '{0}_pipeline.parset'.format(self.name)

        # Define extra parameters needed for this operation
        self.direction.set_imcal_parameters(parset, bands)
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
        loopcount = max(1, self.parset['calibration_specific']['max_selfcal_loops'])

        # Task for smoothing
        if self.parset['calibration_specific']['spline_smooth2d']:
            smooth_amps_task = 'smooth_amps_spline'
        else:
            smooth_amps_task = 'smooth_amps'

        # Parset and sky model for initial solve (not used for later selfcal
        # stages). Note that these must be full paths, as the entries in these
        # steps do not prepend the parset or sky model directories (as is
        # generally done otherwise)
        if self.direction.peel_skymodel is not None:
            initial_selfcal_skymodel = self.direction.peel_skymodel
            initial_selfcal_parset = os.path.join(self.factor_parset_dir,
                'facet_dirdep_phaseonly_solve_skymodel.parset')
        else:
            initial_selfcal_skymodel = os.path.join(self.factor_skymodel_dir,
                'empty.skymodel')
            initial_selfcal_parset = os.path.join(self.factor_parset_dir,
                'facet_dirdep_phaseonly_solve.parset')

        # Parset for slow gain solve
        if self.direction.solve_all_correlations:
            selfcal_caltype = 'fulljones'
            fourpol = True # plot all correlations
        else:
            selfcal_caltype = 'diagonal'
            fourpol = False

        self.parms_dict.update({'ms_files_single': ms_files_single,
                                'ms_files_grouped': str(ms_files),
                                'skymodels': skymodels,
                                'dir_indep_parmDBs': dir_indep_parmDBs,
                                'initial_selfcal_skymodel': initial_selfcal_skymodel,
                                'initial_selfcal_parset': initial_selfcal_parset,
                                'selfcal_caltype': selfcal_caltype,
                                'fourpol': fourpol,
                                'loopcount': loopcount,
                                'smooth_amps_task': smooth_amps_task})

    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object for later use
        self.direction.input_files_single_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'input_files_single.mapfile')
        self.direction.shifted_model_data_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'corrupt_final_model.mapfile')
        self.direction.dir_indep_parmdbs_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'dir_indep_instrument_parmdbs.mapfile')
        self.direction.dir_indep_skymodels_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'full_skymodels.mapfile')
        self.direction.dir_indep_facet_skymodels_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'make_facet_skymodels_all.mapfile')
        self.direction.dir_dep_parmdb_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'merge_selfcal_parmdbs.mapfile')
        self.direction.selfcal_plots_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'make_selfcal_plots.mapfile')
        if self.direction.skip_facet_imaging:
            self.direction.facet_model_mapfile = os.path.join(self.pipeline_mapfile_dir,
                'blank_model.mapfile')
        else:
            self.direction.facet_model_mapfile = os.path.join(self.pipeline_mapfile_dir,
                'final_model_rootnames.mapfile')
        self.direction.facet_image_mapfile = {'facetimage': os.path.join(self.pipeline_mapfile_dir,
            'final_image.mapfile')} # this attribute is a dictionary to allow multiple image mapfiles
        self.direction.facet_premask_mapfile = {'facetimage': os.path.join(self.pipeline_mapfile_dir,
            'premask.mapfile')} # this attribute is a dictionary to allow multiple mask mapfiles
        self.direction.wsclean_modelimg_size_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'pad_model_images.padsize.mapfile')
        self.direction.diff_models_field_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'shift_diff_model_to_field.mapfile')
        self.direction.verify_subtract_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'verify_subtract.break.mapfile')
        self.direction.image_data_mapfile_unconcat = os.path.join(self.pipeline_mapfile_dir,
            'concat_averaged_input.mapfile')
        self.direction.image_data_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'full_image_input.mapfile')
        self.direction.preapply_parmdb_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'create_preapply_parmdb.mapfile')

        # We also need to save the averaging steps for the image_data, so that for
        # any subsequent imaging runs that use these data, we can determine
        # new averaging steps that will give the correct cumulative averaging
        self.direction.full_res_facetimage_freqstep = self.direction.facetimage_freqstep
        self.direction.full_res_facetimage_timestep = self.direction.facetimage_timestep

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

        # Delete all data used only for selfcal as they're no longer needed.
        # Note: we keep the data if selfcal failed verification, so that the user
        # can check them for problems
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'make_sourcedb_all_facet_sources.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'make_sourcedb_cal_facet_sources.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_averaged_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'average_pre_compressed.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'average_post_compressed.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'corrupt_final_model.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'predict_all_model_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'shift_cal.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'final_image1.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'make_concat_corr.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'make_blavg_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'sorted_groups.mapfile_groups'),
            os.path.join(self.pipeline_mapfile_dir, 'average0.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat0_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat1_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat2_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat3_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat4_input.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'wsclean_image01_imagename.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'wsclean_image11_imagename.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'wsclean_image21_imagename.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'wsclean_image31_imagename.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'wsclean_image41_imagename.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'apply_amp1.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'apply_amp2.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'apply_phaseonly1.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'apply_phaseonly2.mapfile')]
        if not self.parset['keep_avg_facet_data'] and self.direction.name != 'target':
            # Add averaged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to reimage by hand (e.g.,
            # with a different weighting). They are always kept for the target
            self.direction.cleanup_mapfiles.extend([
                os.path.join(self.pipeline_mapfile_dir, 'concat_averaged_input.mapfile'),
                os.path.join(self.pipeline_mapfile_dir, 'apply_dir_dep_sorted_groups.mapfile_groups')])
        if not self.parset['keep_unavg_facet_data']:
            # Add unaveraged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to phase shift them to
            # another direction (e.g., to combine several facets together before
            # imaging them all at once)
            self.direction.cleanup_mapfiles.extend([
                os.path.join(self.pipeline_mapfile_dir, 'shift_empty.mapfile'),
                os.path.join(self.pipeline_mapfile_dir, 'sorted_groups_shift_empty.mapfile_groups')])
        if self.direction.selfcal_ok or not self.parset['calibration_specific']['exit_on_selfcal_failure']:
            self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
            self.direction.cleanup()


class FacetPeel(OutlierPeel):
    """
    Operation to peel a facet
    """
    def __init__(self, parset, bands, direction):
        super(FacetPeel, self).__init__(parset, bands, direction,
            name='FacetPeel')

        # Set the pipeline parset to use (the outlierpeel one)
        self.pipeline_parset_template = 'outlierpeel_pipeline.parset'


class FacetSub(Operation):
    """
    Operation to subtract improved model of a facet
    """
    def __init__(self, parset, bands, direction):
        super(FacetSub, self).__init__(parset, bands, direction,
            name='FacetSub')

        # Set the pipeline parset to use
        if self.direction.peel_calibrator:
            # Set parset template to replace parset
            self.pipeline_parset_template = '{0}_replace_pipeline.parset'.format(self.name)
        else:
            self.pipeline_parset_template = '{0}_pipeline.parset'.format(self.name)


    def finalize(self):
        """
        Finalize this operation
        """
        # Delete the full-resolution model-difference files
        if hasattr(self.direction, 'diff_models_field_mapfile'):
            self.direction.cleanup_mapfiles.append(self.direction.diff_models_field_mapfile)
        if hasattr(self.direction, 'subtracted_data_new_mapfile'):
            # Delete the improved subtracted data from peelimage op
            self.direction.cleanup_mapfiles.append(self.direction.subtracted_data_new_mapfile)
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()


class FacetSubReset(Operation):
    """
    Operation to reset the subtraction of improved model of a facet
    """
    def __init__(self, parset, bands, direction):
        super(FacetSubReset, self).__init__(parset, bands, direction,
            name='FacetSubReset')

        # Define extra parameters needed for this operation
        self.direction.set_imcal_parameters(parset, bands)
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
                                'ms_files_grouped': str(ms_files),
                                'skymodels': skymodels,
                                'dir_indep_parmDBs': dir_indep_parmDBs})


    def finalize(self):
        """
        Finalize this operation
        """
        # Delete temp data
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'regroup_shift_empty.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'corrupt_final_model.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'predict_all_model_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'shift_diff_model_to_field.mapfile_groups')]
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()


class FacetImage(Operation):
    """
    Operation to make the full image of a facet

    Two pipelines can be run, depending on whether an improved model (from
    selfcal) exists:

    facetimage_skymodel_pipeline.parset - runs imaging with original
        skymodel

    facetimage_imgmodel_pipeline.parset - runs imaging with improved
        model

    """
    def __init__(self, parset, bands, direction, cellsize_arcsec, robust,
        taper_arcsec, min_uv_lambda):
        # Set name from the facet imaging parameters. If the parameters are all
        # the same as those used for selfcal, just use 'FacetImage'; otherwise,
        # append the parameters
        selfcal_robust = parset['imaging_specific']['selfcal_robust']
        if (cellsize_arcsec != parset['imaging_specific']['selfcal_cellsize_arcsec'] or
            robust != selfcal_robust or taper_arcsec != 0.0 or
            min_uv_lambda != parset['imaging_specific']['selfcal_min_uv_lambda']):
            name = 'FacetImage_c{0}r{1}t{2}u{3}'.format(round(cellsize_arcsec, 1),
                    round(robust, 2), round(taper_arcsec, 1), round(min_uv_lambda, 1))
        else:
            name = 'FacetImage'
        super(FacetImage, self).__init__(parset, bands, direction,
            name=name)

        # Set flag for full-resolution run (used in finalize() to ensure that averaging
        # of the calibrated data is not too much for use by later imaging runs)
        if cellsize_arcsec == parset['imaging_specific']['selfcal_cellsize_arcsec']:
            self.full_res = True
        else:
            self.full_res = False

        # Check whether data files needed for imaging exist already or need to
        # be regenerated
        if hasattr(self.direction, 'image_data_mapfile'):
            self.direction.use_existing_data = self.check_existing_files(self.direction.image_data_mapfile)
            if (self.direction.use_existing_data and
                'facetselfcal' in self.direction.image_data_mapfile and
                self.parset['imaging_specific']['reimage_selfcaled']):
                # Old data exist but are from facetselfcal. Since reimage_selfcal is True,
                # we will not use these data but instead generate new updated ones
                self.log.debug('All files exist but are not up-to-date. They will be regenerated')
                self.direction.use_existing_data = False

                # Store mapfile for old data from selfcal for later clean up
                self.direction.image_data_mapfile_selfcal = self.direction.image_data_mapfile
        else:
            self.direction.use_existing_data = False

        # Set the pipeline parset
        if not self.direction.selfcal_ok:
            # Set parset template to sky-model parset
            self.pipeline_parset_template = 'facetimage_skymodel_pipeline.parset'
        else:
            # Set parset template to facet model-image parset
            self.pipeline_parset_template = 'facetimage_imgmodel_pipeline.parset'

        # Define extra parameters needed for this operation
        self.direction.set_imcal_parameters(parset, bands, cellsize_arcsec, robust,
            taper_arcsec, min_uv_lambda, imaging_only=True,
            use_existing_data=self.direction.use_existing_data,
            existing_data_freqstep=self.direction.full_res_facetimage_freqstep,
            existing_data_timestep=self.direction.full_res_facetimage_timestep)
        if self.direction.use_existing_data:
            self.log.debug('Suitable calibrated data exist and will be used for reimaging')
            # Set flag that determines whether additional averaging is to be done
            if (self.direction.facetimage_freqstep != 1 or self.direction.facetimage_timestep != 1):
                self.direction.average_image_data = True
        else:
            self.direction.image_data_mapfile = 'concat_averaged_compressed_map.output.mapfile'
            self.log.debug('No suitable calibrated data exist for reimaging. They will be generated')

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
        # Add output datamaps to direction object for later use
        self.direction.facet_image_mapfile[self.name.lower()] = os.path.join(self.pipeline_mapfile_dir,
            'final_image.mapfile')
        self.direction.facet_premask_mapfile[self.name.lower()] = os.path.join(self.pipeline_mapfile_dir,
            'premask.mapfile')

        # Store the image_data_mapfile for use by other imaging runs. We do not
        # update this if use_existing_data is True, as in this case it should
        # point to the mapfile of the existing data for this direction. Also, we
        # set this only if this is a full-res imaging run, to ensure that the
        # averaging is not too much for use by later imaging runs
        if not self.direction.use_existing_data and self.full_res:
            self.direction.image_data_mapfile_unconcat = os.path.join(self.pipeline_mapfile_dir,
                'concat_averaged_input.mapfile')
            self.direction.image_data_mapfile = os.path.join(self.pipeline_mapfile_dir,
                'concat_averaged_compressed.mapfile')

            # We also need to save the averaging steps for these data, so that
            # for the next imaging run, we can determine new averaging steps
            # that will give the correct cumulative averaging
            self.direction.full_res_facetimage_freqstep = self.direction.facetimage_freqstep
            self.direction.full_res_facetimage_timestep = self.direction.facetimage_timestep

        # Delete temp data
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'image1.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'corrupt_final_model.mapfile')]
        if ((not self.parset['keep_avg_facet_data'] and self.direction.name != 'target') or
           self.direction.use_existing_data):
            # Add averaged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to reimage by hand (e.g.,
            # with a different weighting) or for subsequent imaging runs. They
            # are always kept for the target direction
            self.direction.cleanup_mapfiles.extend([
                os.path.join(self.pipeline_mapfile_dir, 'concat_averaged_input.mapfile'),
                os.path.join(self.pipeline_mapfile_dir, 'sorted_groups.mapfile_groups')])
        if not self.direction.use_existing_data and hasattr(self.direction, 'image_data_mapfile_selfcal'):
            # Add old data from selfcal to files to be deleted, as we have made new improved versions
            self.direction.cleanup_mapfiles.append(self.direction.image_data_mapfile_selfcal)
        if not self.parset['keep_unavg_facet_data']:
            # Add unaveraged calibrated data for the facet to files to be deleted
            self.direction.cleanup_mapfiles.extend([
                os.path.join(self.pipeline_mapfile_dir, 'shift_empty.mapfile'),
                os.path.join(self.pipeline_mapfile_dir, 'sorted_groups_shift_empty.mapfile_groups')])
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()


class FacetPeelImage(Operation):
    """
    Operation to make the full image of a facet after peeling

    Note, we do not allow cellsize, robust, or taper parameters, as they must
    match the selfcal ones
    """
    def __init__(self, parset, bands, direction):
        super(FacetPeelImage, self).__init__(parset, bands, direction,
            name='FacetPeelImage')

        # Define extra parameters needed for this operation
        self.direction.set_imcal_parameters(parset, bands, imaging_only=True)
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
        # Add output datamaps to direction object for later use
        self.direction.facet_image_mapfile['facetimage'] = os.path.join(self.pipeline_mapfile_dir,
            'final_image.mapfile')
        self.direction.facet_premask_mapfile['facetimage'] = os.path.join(self.pipeline_mapfile_dir,
            'premask.mapfile')
        self.direction.subtracted_data_new_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'subtract_facet_model.mapfile')
        self.direction.facet_model_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'final_model_rootnames.mapfile')
        self.direction.wsclean_modelimg_size_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'pad_model_images.padsize.mapfile')

        # Delete temp data
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'corrupt_final_model.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'concat_averaged_input.mapfile')]
        if not self.parset['keep_avg_facet_data'] and self.direction.name != 'target':
            # Add averaged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to reimage by hand (e.g.,
            # with a different weighting). They are always kept for the target
            self.direction.cleanup_mapfiles.append(
                os.path.join(self.pipeline_mapfile_dir, 'concat_averaged_compressed.mapfile'))
        if not self.parset['keep_unavg_facet_data']:
            # Add unaveraged calibrated data for the facet to files to be deleted
            self.direction.cleanup_mapfiles.append(
                os.path.join(self.pipeline_mapfile_dir, 'shift_empty.mapfile'))
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()
