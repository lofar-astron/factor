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
        if self.direction.contains_target:
            loopcount = max(1, self.parset['calibration_specific']['target_max_selfcal_loops'])
        else:
            loopcount = max(1, self.parset['calibration_specific']['max_selfcal_loops'])

        # Task for smoothing
        if self.parset['calibration_specific']['spline_smooth2d']:
            smooth_amps_task = 'smooth_amps_spline'
        else:
            smooth_amps_task = 'smooth_amps'

        # Set mapfile for selfcal operations
        if self.local_selfcal_scratch_dir is not None:
            if self.direction.pre_average:
                self.direction.phase_concat_data_mapfile = 'make_concat_blavg_data_sync_mapfile.output.mapfile'
            else:
                self.direction.phase_concat_data_mapfile = 'make_concat_data_sync_mapfile.output.mapfile'
            self.direction.concat_data_mapfile = 'make_concat_data_sync_mapfile.output.mapfile'
        else:
            if self.direction.pre_average:
                self.direction.phase_concat_data_mapfile = 'concat_blavg_data.output.mapfile'
            else:
                self.direction.phase_concat_data_mapfile = 'concat_data.output.mapfile'
            self.direction.concat_data_mapfile = 'concat_data.output.mapfile'

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
        self.direction.dir_indep_parmdbs_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'dir_indep_instrument_parmdbs.mapfile')
        self.direction.dir_indep_skymodels_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'full_skymodels.mapfile')
        self.direction.dir_indep_facet_skymodels_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'make_facet_skymodels_all.mapfile')
        self.direction.dir_dep_parmdb_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'merge_selfcal_parmdbs.mapfile')
        self.direction.converted_parmdb_mapfile= os.path.join(self.pipeline_mapfile_dir,
            'convert_merged_selfcal_parmdbs.mapfile')
        self.direction.selfcal_plots_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'make_selfcal_plots.mapfile')
        # Store mapfile of facet image and mask, needed by the fieldmosaic op. Here
        # we store it under the facetimage entry which will get overwritten later
        # if the facet is reimaged
        if hasattr(self.direction, 'facet_image_mapfile'):
            self.direction.facet_image_mapfile['facetimage'] = os.path.join(self.pipeline_mapfile_dir,
                'final_image.mapfile')
            self.direction.facet_premask_mapfile['facetimage'] = os.path.join(self.pipeline_mapfile_dir,
                'premask.mapfile')
        else:
            self.direction.facet_image_mapfile = {'facetimage': os.path.join(self.pipeline_mapfile_dir,
                'final_image.mapfile')}
            self.direction.facet_premask_mapfile = {'facetimage': os.path.join(self.pipeline_mapfile_dir,
                'premask.mapfile')}
        self.direction.wsclean_modelimg_size_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'pad_model_images.padsize.mapfile')
        self.direction.diff_models_field_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'predict_and_difference_models.mapfile')
        self.direction.sourcedb_new_facet_sources = os.path.join(self.pipeline_mapfile_dir,
            'make_sourcedb_new_facet_sources.mapfile')
        self.direction.verify_subtract_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'verify_subtract.break.mapfile')
        if self.direction.create_preapply_parmdb:
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
            os.path.join(self.pipeline_mapfile_dir, 'shift_cal.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'shift_cal_dir_indep.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'make_concat_corr.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'make_blavg_data.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'sorted_groups.mapfile_groups'),
            os.path.join(self.pipeline_mapfile_dir, 'sorted_average0_groups.mapfile_groups'),
            os.path.join(self.pipeline_mapfile_dir, 'average0.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'average2.mapfile'),
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
            os.path.join(self.pipeline_mapfile_dir, 'apply_output.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'apply_phaseonly1.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'apply_phaseonly2.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'subtract_high.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'prepare_imaging_data.mapfile')
            ]
        if self.direction.selfcal_ok or not self.parset['calibration_specific']['exit_on_selfcal_failure']:
            self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
            self.direction.cleanup()
        self.cleanup()


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


    def finalize(self):
        """
        Finalize this operation
        """
        # Delete the full-resolution model-difference files
        if hasattr(self.direction, 'diff_models_field_mapfile'):
            self.direction.cleanup_mapfiles.append(self.direction.diff_models_field_mapfile)
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()
        self.cleanup()


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
        self.direction.cleanup_mapfiles = [os.path.join(self.pipeline_mapfile_dir, 'predict_and_difference_models.mapfile')]
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()
        self.cleanup()


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

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'facetimage_pipeline.parset'

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
        else:
            self.direction.use_existing_data = False

        # Define extra parameters needed for this operation
        self.direction.set_imcal_parameters(parset, bands, cellsize_arcsec, robust,
            taper_arcsec, min_uv_lambda, imaging_only=True,
            use_existing_data=self.direction.use_existing_data,
            existing_data_freqstep=self.direction.full_res_facetimage_freqstep,
            existing_data_timestep=self.direction.full_res_facetimage_timestep)
        if self.direction.use_existing_data:
            self.log.debug('Suitable calibrated data exist and will be used for (re)imaging')
            # Set flag that determines whether additional averaging is to be done
            if (self.direction.facetimage_freqstep != 1 or self.direction.facetimage_timestep != 1):
                self.direction.average_image_data = True
        else:
            self.direction.image_data_mapfile = 'create_compressed_mapfile.output.mapfile'
            self.log.debug('No suitable calibrated data exist for (re)imaging. They will be generated')

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
        # Store mapfile of facet image and mask, needed by the fieldmosaic op
        if hasattr(self.direction, 'facet_image_mapfile'):
            self.direction.facet_image_mapfile[self.name.lower()] = os.path.join(self.pipeline_mapfile_dir,
                'final_image.mapfile')
            self.direction.facet_premask_mapfile[self.name.lower()] = os.path.join(self.pipeline_mapfile_dir,
                'premask.mapfile')
        else:
            self.direction.facet_image_mapfile = {self.name.lower(): os.path.join(self.pipeline_mapfile_dir,
                'final_image.mapfile')}
            self.direction.facet_premask_mapfile = {self.name.lower(): os.path.join(self.pipeline_mapfile_dir,
                'premask.mapfile')}

        # Store the image_data_mapfile for use by other imaging runs. We do not
        # update this if use_existing_data is True, as in this case it should
        # point to the mapfile of the existing data for this direction. Also, we
        # set this only if this is a full-res imaging run, to ensure that the
        # averaging is not too much for use by later imaging runs
        if not self.direction.use_existing_data and self.full_res:
            self.direction.image_data_mapfile = os.path.join(self.pipeline_mapfile_dir,
                'imaging_input.mapfile')

            # We also need to save the averaging steps for these data, so that
            # for the next imaging run, we can determine new averaging steps
            # that will give the correct cumulative averaging
            self.direction.full_res_facetimage_freqstep = self.direction.facetimage_freqstep
            self.direction.full_res_facetimage_timestep = self.direction.facetimage_timestep

        # Delete temp data
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'image1.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'add_all_facet_sources.mapfile'),
            os.path.join(self.pipeline_mapfile_dir, 'corrupt_final_model.mapfile')]
        if ((not self.parset['keep_avg_facet_data'] and self.direction.contains_target) or
           self.direction.use_existing_data):
            # Add averaged calibrated data for the facet to files to be deleted.
            # These are only needed if the user wants to reimage by hand (e.g.,
            # with a different weighting) or for subsequent imaging runs. They
            # are always kept for the target direction
            self.direction.cleanup_mapfiles.extend([
                os.path.join(self.pipeline_mapfile_dir, 'concat_averaged_input.mapfile'),
                os.path.join(self.pipeline_mapfile_dir, 'sorted_groups.mapfile_groups')])
        if not self.parset['keep_unavg_facet_data']:
            # Add unaveraged calibrated data for the facet to files to be deleted
            self.direction.cleanup_mapfiles.extend([
                os.path.join(self.pipeline_mapfile_dir, 'shift_empty.mapfile'),
                os.path.join(self.pipeline_mapfile_dir, 'sorted_groups_shift_empty.mapfile'),
                os.path.join(self.pipeline_mapfile_dir, 'sorted_groups_shift_empty.mapfile_groups')])
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()
        self.cleanup()
