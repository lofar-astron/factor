"""
Module that holds all facet operations

Operations are largely defined by logical groups of actions, but
they also break up the processing into parallel or serial steps.

Classes
-------
FacetAdd : Operation
    Adds facet sources to data
FacetSelfcal : Operation
    Runs the selfcal and imaging of a facet
FacetSub : Operation
    Subtracts all facet sources from data
FacetAddAllFinal : Operation
    Adds all facet sources from final model
FacetImageFinal : Operation
    Images the entire facet

"""
import os
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
        self.parms_dict = {'input_dir': parset['dir_ms'],
                           'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'add_all_parset': add_all_parset,
                           'add_cal_parset': add_cal_parset,
                           'dir_indep_parmdb_name': parset['parmdb_name'],
                           'skymodels': skymodels,
                           'facet_ra': self.direction.ra,
                           'facet_dec': self.direction.dec,
                           'cal_radius_deg': self.direction.cal_radius_deg,
                           'facet_state_file': self.direction.save_file,
                           'hosts': self.node_list,
                           'max_cpus_per_node': self.max_cpus_per_node}


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


class FacetSelfcal(Operation):
    """
    Operation to selfcal one or more directions
    """
    def __init__(self, parset, bands, direction):
        super(FacetSelfcal, self).__init__(parset, bands, direction,
            name='FacetSelfcal')

        # Define parameters needed for this operation
        if self.direction.nchannels > 1:
            nterms = 2
            casa_suffix = '.tt0'
        else:
            nterms = 1
            casa_suffix = None
        self.parms_dict = {'input_dir': parset['dir_ms'],
                           'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'script_dir': self.factor_script_dir,
                           'shifted_cal_bands_datamap': self.direction.shifted_cal_bands_datamap,
                           'shifted_all_bands_datamap': self.direction.shifted_all_bands_datamap,
                           'shifted_empty_bands_datamap': self.direction.shifted_empty_bands_datamap,
                           'dir_indep_parmdbs_datamap': self.direction.dir_indep_parmdbs_datamap,
                           'field_ra': self.direction.field_ra,
                           'field_dec': self.direction.field_dec,
                           'wplanes': self.direction.wplanes,
                           'casa_suffix': casa_suffix,
                           'facet_imsize': self.direction.facet_imsize,
                           'cal_wplanes': self.direction.wplanes,
                           'cal_imsize': self.direction.cal_imsize,
                           'nterms': nterms,
                           'nchannels': self.direction.nchannels,
                           'chunk_width': (self.direction.solint_a-1)*2,
                           'solint_p': self.direction.solint_p,
                           'solint_a': self.direction.solint_a,
                           'facet_state_file': self.direction.save_file,
                           'hosts': self.direction.hosts,
                           'max_cpus_per_node': self.max_cpus_per_node}


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamap to direction object
        self.direction.shifted_cal_concat_datamap = os.path.join(self.mapfile_dir,
            'shifted_cal_concat.datamap')
        self.direction.dir_dep_parmdb_datamap = os.path.join(self.mapfile_dir,
            'dir_dep_parmdb.datamap')
        self.direction.facet_image_mapfile = os.path.join(self.mapfile_dir,
            'facet_image.datamap')
        try:
            ok_datamap = DataMap.load(os.path.join(self.mapfile_dir,
                'verify_subtract.datamap'))
            self.direction.selfcal_ok = ok_datamap[0].item
        except:
            pass


class FacetSub(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, bands, direction):
        super(FacetSub, self).__init__(parset, bands, direction,
            name='FacetSub')

        # Define parameters needed for this operation
        self.parms_dict = {'input_dir': parset['dir_ms'],
                           'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'shifted_all_bands_datamap': self.direction.shifted_all_bands_datamap,
                           'dir_dep_parmdb_datamap': self.direction.dir_dep_parmdb_datamap,
                           'input_bands_datamap': self.direction.input_bands_datamap,
                           'field_ra': self.direction.field_ra,
                           'field_dec': self.direction.field_dec,
                           'hosts': self.direction.hosts}


    def finalize(self):
        """
        Finalize this operation
        """
        pass


class FacetAddAllFinal(Operation):
    """
    Operation to add all sources in the facet to data (final)
    """
    def __init__(self, parset, bands, direction):
        super(FacetAddAllFinal, self).__init__(parset, bands, direction,
            name=name)


def FacetImageFinal(FacetImage):
    """
    Operation to make final facet image
    """
    def __init__(self, parset, bands, direction):
        super(FacetImageFinal, self).__init__(parset, bands, direction,
            name='FacetImageFinal')

