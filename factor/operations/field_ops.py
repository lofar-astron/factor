"""
Module that holds all field (non-facet-specific) operations

Classes
-------
InitSubtract : Operation
    Images each band at high and low resolution to make and subtract sky models
MakeMosaic : Operation
    Makes a mosaic of the field from the facet images

"""
import os
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap

class InitSubtract(Operation):
    """
    Operation to create empty datasets
    """
    def __init__(self, parset, bands, direction):
        super(InitSubtract, self).__init__(parset, bands, direction,
            name='InitSubtract')

        # Define parameters needed for this operation
        self.parms_dict = {'input_dir': parset['dir_ms'],
                           'npix_high': self.direction.imsize_high_res,
                           'npix_low': self.direction.imsize_low_res,
                           'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'dir_indep_parmdb_name': parset['parmdb_name'],
                           'hosts': self.node_list,
                           'max_cpus_per_node': self.max_cpus_per_node}


    def finalize(self):
        """
        Finalize this operation
        """
        # Add skymodels to band objects
        merged_skymodel_datamap = os.path.join(self.mapfile_dir,
            'merged_skymodels.datamap')
        datamap = DataMap.load(merged_skymodel_datamap)
        for band, item in zip(self.bands, datamap):
            band.skymodel_dirindep = item.file
        self.cleanup_mapfiles.append(os.path.join(self.mapfile_dir,
            'averaged_data.datamap'))


class MakeMosaic(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, bands, direction):
        super(MakeMosaic, self).__init__(parset, bands, direction,
            name='MakeMosaic')

        # Define parameters needed for this operation
        self.parms_dict = {'parset_dir': self.factor_parset_dir,
                           'mapfile_dir': self.mapfile_dir}
