"""
Module that holds all field (non-facet-specific) operations

Classes
-------
InitSubtract : Operation
    Images each band at high and low resolution to make and subtract sky models
MakeMosaic : Operation
    Makes a mosaic from the facet images

"""
import os
from factor.lib.operation import Operation


class InitSubtract(Operation):
    """
    Operation to create empty datasets
    """
    def __init__(self, parset, bands, direction):
        super(InitSubtract, self).__init__(parset, bands, direction,
            name='InitSubtract')

        # Define parameters needed for this operation
        self.parms_dict = {'input_dir': parset['dir_ms'],
                           'npix_high': direction.imsize_high_res,
                           'npix_low': direction.imsize_low_res,
                           'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'dir_indep_parmdb_name': parset['parmdb_name']}

        # Add info to direction object
        merged_skymodel_datamap = os.path.join(self.mapfile_dir,
            'merged_skymodels.datamap')
        self.direction.merged_skymodel_datamap = merged_skymodel_datamap


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
