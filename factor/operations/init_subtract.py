"""
Operation: InitSubtract

Images each band at high and low resolution to make and subtract sky models.
"""
import os
from factor.lib.operation import Operation
from factor.actions.images import MakeImage
from factor.actions.models import MakeSkymodelFromModelImage, MergeSkymodels
from factor.actions.calibrations import Subtract

class InitSubtract(Operation):
    """
    Operation to create empty datasets
    """"
    def __init__(self, parset, bands, direction=None, reset=False):
        super(InitSubtract, self).__init__(parset, bands, direction=direction,
            reset=reset, name='InitSubtract')


    def run_steps(self):
        """
        Run the steps for this operations
        """
        from factor.operations.hardcoded_param import init_subtract_test_quick as p

        bands = self.bands

        if 'dir_node' in self.parset:
            localdir = self.parset['dir_node']
            self.log.info('Using {0} for imaging...'.format(localdir))
        else:
            localdir = None

        # Check operation state
        if os.path.exists(self.statebasename+'.done'):
            for band in bands:
                band.skymodel_dirindep = ''
            self.finalize()
            return

        # Make initial data map for the empty datasets
        subtracted_all_datamap = self.make_datamap([band.file for band in bands],
            'initial_subtracted')

        self.log.info('High-res imaging...')
        action = MakeImage(self.parset, subtracted_all_datamap,
            p['imagerh'], prefix='init_highres', localdir=localdir)
        high_res_images_datamap, high_res_models_mapfile, high_res_masks_mapfile = \
            action.get_results()
