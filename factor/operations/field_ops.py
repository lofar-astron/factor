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
    def __init__(self, parset, bands, reset=False):
        super(InitSubtract, self).__init__(parset, bands, direction=None,
            reset=reset, name='InitSubtract')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.images import MakeImage
        from factor.actions.models import MakeSkymodelFromModelImage, MergeSkymodels
        from factor.actions.calibrations import Subtract
        from factor.actions.visibilities import Average
        from factor.lib.datamap_lib import write_mapfile, read_mapfile
        from factor.operations.hardcoded_param import init_subtract_test_quick as p

        bands = self.bands

        if 'dir_node' in self.parset:
            localdir = self.parset['dir_node']
            self.log.info('Using {0} for imaging...'.format(localdir))
        else:
            localdir = None

        # Check operation state
        if os.path.exists(self.statebasename+'.done'):
            merged_skymodels_mapfile = '' # fill in!!
            for band, skymodel in zip(bands, read_mapfile(merged_skymodel_mapfile)):
                band.skymodel_dirindep = skymodel
            return

        # Make initial data maps for the empty datasets and their dir-indep
        # instrument parmdbs
        subtracted_all_mapfile = write_mapfile([band.file for band in bands],
            self.name, prefix='subtracted_all', working_dir=self.parset['dir_working'])
        dir_indep_parmdbs_mapfile = write_mapfile([os.path.join(band.file,
            band.dirindparmdb) for band in bands], self.name,
            prefix='dir_indep_parmdbs', working_dir=self.parset['dir_working'])

        self.log.info('High-res imaging...')
        action = MakeImage(self.parset, subtracted_all_mapfile, p['imagerh'],
            prefix='highres_image', localdir=localdir)
        highres_image_basenames_mapfile = action.run()

        self.log.info('Making high-res sky model...')
        action = MakeSkymodelFromModelImage(self.parset, highres_image_basenames_mapfile,
            p['modelh'], prefix='highres_model')
        highres_skymodels_mapfile = action.run()

        self.log.info('Subtracting high-res sky model...')
        action = Subtract(self.parset, subtracted_all_mapfile, p['calibh'],
            model_datamap=highres_skymodels_mapfile,
            parmdb_datamap=dir_indep_parmdbs_mapfile, prefix='highres_subtract')
        action.run()

        self.log.info('Averaging...')
        action = Average(self.parset, subtracted_all_mapfile, p['avgl'],
            prefix='highres_average')
        avg_files_mapfile = action.run()

        self.log.info('Low-res imaging...')
        action = MakeImage(self.parset, avg_files_mapfile, p['imagerl'],
            prefix='lowres_image', localdir=localdir)
        lowres_image_basenames_mapfile = action.run()

        self.log.info('Making low-res sky model...')
        action = MakeSkymodelFromModelImage(self.parset, lowres_image_basenames_mapfile,
            p['modell'], prefix='lowres_model')
        lowres_skymodels_mapfile = action.run()

        self.log.info('Subtracting low-res sky model...')
        action = Subtract(self.parset, subtracted_all_mapfile, p['calibl'],
            model_datamap=lowres_skymodels_mapfile,
            parmdb_datamap=dir_indep_parmdbs_mapfile, prefix='lowres_subtract')
        action.run()

        self.log.info('Merging low- and high-res sky models...')
        action = MergeSkymodels(self.parset, [low_res_skymodels_mapfile,
            high_res_skymodels_mapfile], p['merge'], prefix='init_merge')
        merged_skymodels_mapfile = action.run()
        for band, skymodel in zip(bands, read_mapfile(merged_skymodel_mapfile)):
            band.skymodel_dirindep = skymodel


class MakeMosaic(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(MakeMosaic, self).__init__(parset, bands, direction=direction,
            reset=reset, name='MakeMosaic')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        pass
