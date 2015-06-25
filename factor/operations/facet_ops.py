"""
Module that holds all facet operations

The facet calibration steps are divided into operations on the basis of whether
or not they can be run in parallel or in series.

Classes
-------
FacetAdd : Operation
    Adds facet sources to data. Must be run in series as writes are made to
    original datasets
FacetSelfcal : Operation
    Runs the selfcal and imaging of a facet. May be run in parallel
FacetSub : Operation
    Subtracts all facet sources from data. Must be run in series as writes are
    made to original datasets
FacetAddFinal : Operation
    Adds all facet sources from final model. Must be run in series as writes are
    made to original datasets
FacetImageFinal : Operation
    Images the entire facet. May be run in parallel

"""
import os
import ast
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap


class FacetAdd(Operation):
    """
    Operation to add calibrator source to data
    """
    def __init__(self, parset, bands, direction):
        super(FacetAdd, self).__init__(parset, bands, direction,
            name='FacetAdd')

        # Define parameters needed for this operation
        skymodels = [band.skymodel_dirindep for band in self.bands]
        if self.direction.use_new_sub_data:
            add_all_parset = 'facet_dirindep_add_all_new.parset'
            add_cal_parset = 'facet_dirindep_add_cal_new.parset'
        else:
            add_all_parset = 'facet_dirindep_add_all.parset'
            add_cal_parset = 'facet_dirindep_add_cal.parset'
        self.parms_dict.update({'input_dir': parset['dir_ms'],
                                'add_all_parset': add_all_parset,
                                'add_cal_parset': add_cal_parset,
                                'dir_indep_parmdb_name': parset['parmdb_name'],
                                'skymodels': skymodels,
                                'facet_ra': self.direction.ra,
                                'facet_dec': self.direction.dec,
                                'cal_radius_deg': self.direction.cal_size_deg/2.0,
                                'facet_state_file': self.direction.save_file})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object
        self.direction.input_bands_datamap = os.path.join(self.mapfile_dir,
            'input_bands.datamap')
        self.direction.shifted_all_bands_datamap = os.path.join(self.mapfile_dir,
            'shifted_all_bands.datamap')
        self.direction.shifted_cal_bands_datamap = os.path.join(self.mapfile_dir,
            'shifted_cal_bands.datamap')
        self.direction.shifted_empty_bands_datamap = os.path.join(self.mapfile_dir,
            'shifted_empty_bands.datamap')
        self.direction.dir_indep_parmdbs_datamap = os.path.join(self.mapfile_dir,
            'dir_indep_instrument_parmdbs.datamap')
        self.direction.cleanup_mapfiles.extend([self.direction.shifted_all_bands_datamap,
            self.direction.shifted_cal_bands_datamap,
            self.direction.shifted_empty_bands_datamap])


class FacetSelfcal(Operation):
    """
    Operation to selfcal one or more directions
    """
    def __init__(self, parset, direction):
        super(FacetSelfcal, self).__init__(parset, None, direction,
            name='FacetSelfcal')

        # Define parameters needed for this operation
        if self.direction.nchannels > 1:
            nterms = 2
            casa_suffix = '.tt0'
            wsclean_suffix = '-MFS-image.fits'
        else:
            nterms = 1
            casa_suffix = None
            wsclean_suffix = '-image.fits'
        self.parms_dict.update({'shifted_cal_bands_datamap': self.direction.shifted_cal_bands_datamap,
                                'shifted_all_bands_datamap': self.direction.shifted_all_bands_datamap,
                                'shifted_empty_bands_datamap': self.direction.shifted_empty_bands_datamap,
                                'dir_indep_parmdbs_datamap': self.direction.dir_indep_parmdbs_datamap,
                                'field_ra': self.direction.field_ra,
                                'field_dec': self.direction.field_dec,
                                'wplanes': self.direction.wplanes,
                                'casa_suffix': casa_suffix,
                                'wsclean_suffix': wsclean_suffix,
                                'facet_imsize': self.direction.facet_imsize,
                                'cal_wplanes': self.direction.wplanes,
                                'cal_imsize': self.direction.cal_imsize,
                                'nterms': nterms,
                                'nchannels': self.direction.nchannels,
                                'chunk_width': (self.direction.solint_a-1)*2,
                                'solint_p': self.direction.solint_p,
                                'solint_a': self.direction.solint_a,
                                'region_selfcal': self.direction.region_selfcal,
                                'facet_state_file': self.direction.save_file,
                                'hosts': self.direction.hosts})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamap to direction object
        self.direction.dir_dep_parmdb_datamap = os.path.join(self.mapfile_dir,
            'dir_dep_parmdb.datamap')
        self.direction.facet_image_mapfile = os.path.join(self.mapfile_dir,
            'final_image.datamap')
        self.direction.facet_model_mapfile = os.path.join(self.mapfile_dir,
            'final_model_rootnames.datamap')
        self.verify_subtract_OK_mapfile = os.path.join(self.mapfile_dir,
            'verify_subtract_OK.datamap')
        self.direction.cleanup_mapfiles.extend([os.path.join(self.mapfile_dir,
            'chunk_files.datamap'), os.path.join(self.mapfile_dir,
            'concat1_input.datamap'), os.path.join(self.mapfile_dir,
            'concat2_input.datamap'), os.path.join(self.mapfile_dir,
            'concat3_input.datamap'), os.path.join(self.mapfile_dir,
            'concat4_input.datamap')])

        # Store results of verify_subtract check
        if os.path.exists(self.verify_subtract_OK_mapfile):
            ok_datamap = DataMap.load(self.verify_subtract_OK_mapfile)
            ok_flags = [ast.literal_eval(item.file) for item in ok_datamap]
            if self.parset['strict_selfcal_verification']:
                if all(ok_flags):
                    self.direction.selfcal_ok = True
                else:
                    self.direction.selfcal_ok = False
            else:
                if any(ok_flags):
                    self.direction.selfcal_ok = True
                else:
                    self.direction.selfcal_ok = False


class FacetSub(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, direction):
        super(FacetSub, self).__init__(parset, None, direction,
            name='FacetSub')

        # Define parameters needed for this operation
        self.parms_dict.update({'shifted_all_bands_datamap': self.direction.shifted_all_bands_datamap,
                                'dir_dep_parmdb_datamap': self.direction.dir_dep_parmdb_datamap,
                                'input_bands_datamap': self.direction.input_bands_datamap,
                                'field_ra': self.direction.field_ra,
                                'field_dec': self.direction.field_dec})


    def finalize(self):
        """
        Finalize this operation
        """
        self.direction.facet_model_data_mapfile = os.path.join(self.mapfile_dir,
            'shifted_models.datamap')


class FacetAddFinal(Operation):
    """
    Operation to add all sources in the facet in preparation for final imaging
    """
    def __init__(self, parset, bands, direction):
        super(FacetAddFinal, self).__init__(parset, bands, direction,
            name='FacetAddFinal')

        # Define parameters needed for this operation
        if not self.direction.selfcal_ok:
            # Use dir-indep CC sky models
            skymodels = [band.skymodel_dirindep for band in self.bands]

            # Set parset template to sky-model parset
            self.pipeline_parset_template = '{0}_cc_skymodel_pipeline.parset'.format(self.name)

            self.parms_dict.update({'input_dir': parset['dir_ms'],
                                'dir_indep_parmdb_name': parset['parmdb_name'],
                                'skymodels': skymodels,
                                'facet_ra': self.direction.ra,
                                'facet_dec': self.direction.dec,
                                'facet_state_file': self.direction.save_file})
        else:
            # Set parset template to facet model-image parset
            self.pipeline_parset_template = '{0}_model_image_pipeline.parset'.format(self.name)

            self.parms_dict.update({'input_bands_datamap': self.direction.input_bands_datamap,
                                'dir_dep_parmdb_mapfile': self.direction.dir_dep_parmdb_datamap,
                                'facet_model_data_mapfile': self.direction.facet_model_data_mapfile,
                                'facet_ra': self.direction.ra,
                                'facet_dec': self.direction.dec,
                                'facet_state_file': self.direction.save_file})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object
        self.direction.shifted_all_final_bands_datamap = os.path.join(self.mapfile_dir,
            'shifted_all_final_bands.datamap')
        self.direction.cleanup_mapfiles.extend([self.direction.shifted_all_final_bands_datamap])



class FacetImageFinal(Operation):
    """
    Operation to make final facet image
    """
    def __init__(self, parset, direction):
        super(FacetImageFinal, self).__init__(parset, None, direction,
            name='FacetImageFinal')

        # Define parameters needed for this operation
        if self.direction.nchannels > 1:
            wsclean_suffix = '-MFS-image.fits'
        else:
            wsclean_suffix = '-image.fits'
        self.parms_dict.update({'dir_dep_parmdb_datamap': self.direction.dir_dep_parmdb_datamap,
                                'shifted_all_final_bands_datamap': self.direction.shifted_all_final_bands_datamap,
                                'wsclean_suffix': wsclean_suffix,
                                'facet_imsize': self.direction.facet_imsize,
                                'nchannels': self.direction.nchannels,
                                'facet_state_file': self.direction.save_file,
                                'hosts': self.direction.hosts})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object
        self.direction.facet_image_mapfile = os.path.join(self.mapfile_dir,
            'final_image.datamap')
        self.direction.facet_model_mapfile = os.path.join(self.mapfile_dir,
            'final_model_rootnames.datamap')
