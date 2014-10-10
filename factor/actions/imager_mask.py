"""
Action: imager_mask
Make an image, then extract a mask and re-make the image
Used parameters:
* those of IMAGER action
* those of MAKE_MASK action
Return:
* image filename (in a vector to deal with nterm>1)
* model filename (in a vector to deal with nterm>1)
* mask filename
"""

from factor.lib.action import action

class imager_mask(action):
    """
    Implment the imager_mask action
    """

    def __init__(self, op_name, ms, p, prefix = None, direction = None, clean=True):
        super(imager_mask, self).__init__(op_name, name = 'imager_mask')
        self.ms = ms
        self.p = p
        self.prefix = prefix
        self.direction = direction
        self.clean = clean

    def run(self):

        # this action calls other actions
        import factor.actions as a
        imager = a.imager.imager(self.op_name, self.ms, self.p, self.prefix, self.direction, clean=self.clean)
        imager.run()
        p['image'], _ = imager.get_results()

        make_mask = a.make_mask.make_mask(self.op_name, self.ms, self.p, self.prefix, self.direction, clean=self.clean)
        make_mask.run()
        self.mask = make_mask.get_results()
        p['mask'] = self.mask

        imager.mask = a.imager.imager(self.op_name, self.ms, self.p, self.prefix, self.direction, clean=self.clean)
        imager.run()
        self.image, self.model = imager.get_results()

    def get_results(self):
        """
        Return image, model and mask name
        """
        return self.image, self.model, self.mask

