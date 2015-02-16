"""
Operation: example
An example operation, to be used as a template
"""

from factor.lib.operation import Operation
from factor.lib.datamap import make_datamap
from factor.actions.images import make_image

class example_operation(Operation):

    def run(self):

        from factor.operations.hardcoded_param import example_parms as p

        # Make data map from input bands
        bands = self.bands
        subtracted_datamap = make_datamap(bands)

        # this is a typical pipeline call
        self.log.info('Starting example procedure...')
        image_datamap = make_image(subtracted_datamap, p['image'])

