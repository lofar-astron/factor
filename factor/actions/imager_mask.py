"""
Action: imager_mask
Make an image, then extract a mask and re-make the image
"""

from factor.lib.action import action

class imager_mask(action):
    """
    Implment the imager_mask action
    """

    def __init__(self, op_name, ms, prefix='', niter=None, imsize=None, cell=None, uvrange=None, clean=True):
        super(phase_shifter, self).__init__(op_name, name = 'imager_mask')
        self.ms = ms
        self.prefix = prefix
        self.niter = niter
        self.imsize = imsize
        self.cell = cell
        self.uvrange = uvrange
        self.clean = clean
        self.imagename = imagename = makeimagename(ms, prefix)

    def run(self):
        # this action calls other actions
        import factor.actions as a

        imager = a.imager(self.op_name, self.ms, prefix=self.prefix, niter=self.niter, imsize=self.imsize, \
                            cell=self.cell, uvrange=self.uvrange, mask=None, clean=self.clean)
        imager.run()

        make_mask = a.make_mask(ms)
        make_mask.run()
        self.mask = make_mask.get_results()

        imager.mask = self.mask
        imager.run()
        self.image, self.model = imager.get_results()

    def get_results(self):
        """
        Return image, model and mask name
        """
        return self.image, self.model, self.mask

