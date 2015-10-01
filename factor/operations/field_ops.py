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

    Two pipelines can be run, depending on whether the skymodels are present for
    all bands or not:

    initsubtract_pipeline.parset - runs the full initial subtraction when one
        or more bands lack skymodels

    initsubtract_subonly_pipeline.parset - runs only a subtract step when all
        bands have skymodels

    """
    def __init__(self, parset, bands, direction):
        super(InitSubtract, self).__init__(parset, bands, direction,
            name='InitSubtract')

        # Set the pipeline parset to use
        if all([b.skymodel_dirindep is not None for b in self.bands]):
            self.pipeline_parset_template = '{0}_subonly_pipeline.parset'.format(self.name)
        else:
            self.pipeline_parset_template = '{0}_pipeline.parset'.format(self.name)

        # Define extra parameters needed for this operation (beyond those
        # defined in the master Operation class and as attributes of the
        # direction object)
        input_bands = [b.file for b in self.bands]
        highres_image_sizes = ['{0} {0}'.format(b.imsize_high_res) for b in self.bands]
        lowres_image_sizes = ['{0} {0}'.format(b.imsize_low_res) for b in self.bands]
        skymodels = [band.skymodel_dirindep for band in self.bands]
        dir_indep_parmdbs = [band.dirindparmdb for band in self.bands]
        self.parms_dict.update({'input_bands': input_bands,
                                'highres_image_sizes' : highres_image_sizes,
                                'lowres_image_sizes' : lowres_image_sizes,
                                'skymodels': skymodels,
                                'dir_indep_parmdbs': dir_indep_parmdbs})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add skymodels to band objects if any lack them
        if any([b.skymodel_dirindep is None for b in self.bands]):
            merged_skymodel_datamap = os.path.join(self.mapfile_dir,
                'merged_skymodels.datamap')
            if os.path.exists(merged_skymodel_datamap):
                datamap = DataMap.load(merged_skymodel_datamap)
                for band, item in zip(self.bands, datamap):
                    band.skymodel_dirindep = item.file
                    band.skip = item.skip
            else:
                for band self.bands:
                    band.skymodel_dirindep = None

        # Delete averaged data as they're no longer needed
        self.direction.cleanup_mapfiles = [os.path.join(self.mapfile_dir,
            'averaged_data.datamap')]
        self.direction.cleanup()


class MakeMosaic(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, direction):
        super(MakeMosaic, self).__init__(parset, None, direction,
            name='MakeMosaic')

        self.parms_dict.update({'input_dir': parset['dir_ms']})
