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

        # specify  the image parameters here
        cellsize_highres_deg = 0.000834
        cellsize_lowres_deg = 0.00694
        fieldsize_highres = 2.5
        fieldsize_lowres = 6.5
        maxlambda_highres = 7000
        maxlambda_lowres = 2000

        # Define extra parameters needed for this operation (beyond those
        # defined in the master Operation class and as attributes of the
        # direction object)
        input_files = [b.files for b in self.bands]
        input_files_single = []
        for bandfiles in input_files:
            for filename in bandfiles:
                input_files_single.append(filename)
        dir_indep_parmDBs = []
        for band in self.bands:
            band.set_image_sizes(cellsize_highres_deg=cellsize_highres_deg,cellsize_lowres_deg=cellsize_lowres_deg,
                                 fieldsize_highres=fieldsize_highres,fieldsize_lowres=fieldsize_lowres)
            for parmdb in band.dirindparmdbs:
                dir_indep_parmDBs.append(parmdb)
        band_names = [b.name for b in self.bands]
        highres_image_sizes = ['{0} {0}'.format(b.imsize_high_res) for b in self.bands]
        lowres_image_sizes = ['{0} {0}'.format(b.imsize_low_res) for b in self.bands]
        #skymodels = [band.skymodel_dirindep for band in self.bands]
        self.parms_dict.update({'input_files_single': input_files_single,
                                'input_files_grouped' : str(input_files),
                                'highres_image_sizes' : highres_image_sizes,
                                'lowres_image_sizes' : lowres_image_sizes,
                                'cellsize_highres_deg' : cellsize_highres_deg,
                                'cellsize_lowres_deg' : cellsize_lowres_deg,
                                'maxlambda_highres' : maxlambda_highres,
                                'maxlambda_lowres' : maxlambda_lowres,
                                #'skymodels': skymodels,
                                'dir_indep_parmDBs': dir_indep_parmDBs})


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
                # this should continue to work, but I'm not sure
                assert len(self.bands) == len(datamap)
                for band, item in zip(self.bands, datamap):
                    band.skymodel_dirindep = item.file
                    band.skip = item.skip                
            else:
                for band in self.bands:
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
