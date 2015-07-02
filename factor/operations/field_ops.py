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
        input_bands = [b.file for b in self.bands]
        highres_image_sizes = ['{0} {0}'.format(b.imsize_high_res) for b in self.bands]
        lowres_image_sizes = ['{0} {0}'.format(b.imsize_low_res) for b in self.bands]
        self.parms_dict.update({'input_bands': input_bands,
                                'highres_image_sizes' : highres_image_sizes,
                                'lowres_image_sizes' : lowres_image_sizes,
                                'max_percent_memory' : self.max_percent_memory,
                                'dir_indep_parmdb_name': parset['parmdb_name']})


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
        self.direction.cleanup_mapfiles.append(os.path.join(self.mapfile_dir,
            'averaged_data.datamap'))


class MakeMosaic(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, direction):
        super(MakeMosaic, self).__init__(parset, None, direction,
            name='MakeMosaic')

        # Define parameters needed for this operation
        self.parms_dict.update({'facet_image_filenames': self.direction.facet_image_filenames,
                                'facet_vertices_filenames': self.direction.facet_vertices_filenames})
