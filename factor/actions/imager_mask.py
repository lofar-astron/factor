"""
Action: imager_mask
Make an image, then extract a mask and re-make the image
"""

from factor.lib.action import action

class imager_mask(action):
    """
    Implment the imager_mask action
    """

    def __init__(self, op_name, ms, p, clean=True):
        super(imager_mask, self).__init__(op_name, name = 'imager_mask')
        self.ms = ms
        self.p = p
        self.clean = clean

    def run(self):
        import factor.actions as a
        
        # this action calls other actions
        imager = a.imager.imager(self.op_name, self.ms, self.p, clean=self.clean)
        imager.run()

        make_mask = a.make_mask.make_mask(self.ms)
        make_mask.run()
        self.mask = make_mask.get_results()
        p['mask'] = self.mask

        imager.mask = a.imager.imager(self.op_name, self.ms, self.p, clean=self.clean)
        imager.run()
        self.image, self.model = imager.get_results()

    def get_results(self):
        """
        Return image, model and mask name
        """
        return self.image, self.model, self.mask

