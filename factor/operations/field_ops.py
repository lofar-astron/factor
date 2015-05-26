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
from factor.lib.scheduler_mp import Scheduler


class InitSubtract(Operation):
    """
    Operation to create empty datasets
    """
    def __init__(self, parset, bands, reset=False):
        super(InitSubtract, self).__init__(parset, bands, direction=None,
            reset=reset, name='InitSubtract')

        # Set up imager scheduler (runs at most num_nodes imagers in parallel)
        num_nodes = len(self.parset['cluster_specific']['node_list'])
        self.s = Scheduler(max_procs=num_nodes, name=self.name,
            op_parset=self.parset)


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.images import image_with_mask
        from factor.actions.models import MakeSkymodelFromModelImage, MergeSkymodels, FFT
        from factor.actions.calibrations import Subtract
        from factor.actions.visibilities import Average, ChgCentre
        from factor.lib.datamap_lib import read_mapfile
        from factor.operations.hardcoded_param import init_subtract as p
        from factor.lib.operation_lib import make_chunks, merge_chunks
        import numpy as np

        bands = self.bands

        # Check state
        if self.check_completed(bands):
            return

        if os.path.exists(self.statebasename+'.done'):
            merged_skymodels_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/InitSubtract/MergeSkymodels/merge_output.datamap')
            skymodels, _ = read_mapfile(merged_skymodels_mapfile)
            for band, skymodel in zip(bands, skymodels):
                band.skymodel_dirindep = skymodel
            input_data_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/InitSubtract/input_data.datamap')
            files, _ = read_mapfile(input_data_mapfile)
            for band, f in zip(bands, files):
                band.file = f
            return

        # Make initial data maps for the input datasets and their dir-indep
        # instrument parmdbs.
        input_data_mapfile = self.write_mapfile([band.file for band in bands],
        	prefix='input_data')
        dir_indep_parmdbs_mapfile = self.write_mapfile([band.dirindparmdb for band
        	in bands], prefix='dir_indep_parmdbs')

        self.log.info('High-res imaging...')
        if self.parset['use_chgcentre']:
            self.log.debug('Changing center to zenith...')
            action = ChgCentre(self.parset, input_data_mapfile, {},
                prefix='highres', band=band)
            chgcentre_data_mapfile = self.s.run(action)
            input_to_imager_mapfile = chgcentre_data_mapfile
        else:
            input_to_imager_mapfile = input_data_mapfile
        highres_image_basenames_mapfile = image_with_mask(self, p['imagerh'],
            'highres', [input_to_imager_mapfile])

        self.log.info('Making high-res sky model...')
        action = MakeSkymodelFromModelImage(self.parset,
            highres_image_basenames_mapfile, p['imagerh'], prefix='highres')
        highres_skymodels_mapfile = self.s.run(action)

        self.log.debug('FFTing high-res model image...')
        action = FFT(self.parset, input_data_mapfile,
            highres_image_basenames_mapfile, p['imagerh'], prefix='highres')
        self.s_imager.run(action)

        self.log.info('Subtracting high-res sky model...')
        if self.parset['use_ftw']:
            highres_skymodels_mapfile = None
        action = Subtract(self.parset, input_data_mapfile, p['calibh'],
            None, dir_indep_parmdbs_mapfile, prefix='highres', band=band)
        self.s.run(action)

        self.log.info('Averaging...')
        if self.parset['use_chgcentre']:
            self.log.debug('Changing center to zenith...')
            action = ChgCentre(self.parset, input_data_mapfile, {},
                prefix='lowres')
            chgcentre_data_mapfile = self.s.run(action)
            input_to_avg_mapfile = chgcentre_data_mapfile
        else:
            input_to_avg_mapfile = input_data_mapfile
        action = Average(self.parset, input_to_avg_mapfile, p['avgl'],
            prefix='highres')
        avg_files_mapfile = self.s.run(action)

        self.log.info('Low-res imaging...')
        lowres_image_basenames_mapfile = image_with_mask(self, p['imagerl'],
            'lowres', [avg_files_mapfile])

        self.log.info('Making low-res sky model...')
        action = MakeSkymodelFromModelImage(self.parset, lowres_image_basenames_mapfile,
            p['imagerl'], prefix='lowres')
        lowres_skymodels_mapfile = self.s.run(action)

        self.log.debug('FFTing low-res model image...')
        action = FFT(self.parset, input_data_mapfile, lowres_image_basenames_mapfile,
            p['imagerl'], prefix='lowres')
        self.s_imager.run(action)

        self.log.info('Subtracting low-res sky model...')
        action = Subtract(self.parset, input_data_mapfile, p['calibl'], None,
            dir_indep_parmdbs_mapfile, prefix='lowres')
        self.s.run(action)

        self.log.info('Merging low- and high-res sky models...')
        action = MergeSkymodels(self.parset, lowres_skymodels_mapfile,
            highres_skymodels_mapfile, p['merge'], prefix='merge')
        merged_skymodels_mapfile = self.s.run(action)
        skymodels, _ = read_mapfile(merged_skymodels_mapfile)
        for band, skymodel in zip(bands, skymodels):
            band.skymodel_dirindep = skymodel

        # Save state
        self.set_completed(bands)


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
