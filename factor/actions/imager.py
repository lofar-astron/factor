"""
Action: imager
Make an image
"""

import os
from factor.lib.action import action
from factor.lib.action_lib import makeimagename

class imager(action):
    """
    Implment the imager action
    """

    def __init__(self, op_name, ms, prefix='', niter=None, imsize=None, cell=None, uvrange=None, mask=None, clean=True):
        super(phase_shifter, self).__init__(op_name, name = 'imager')
        self.ms = ms
        self.prefix = prefix
        self.niter = niter
        self.imsize = imsize
        self.cell = cell
        self.uvrange = uvrange
        self.mask = mask
        self.clean = clean
        self.image = makeimagename(ms, prefix)

    def run(self):
        # TODO: implement the template
        template_imager = make_template(ms = self.ms, image = self.image, \
                                niter = self.niter, imsize = self.imsize, cell = self.cell, uvrange = self.uvrange, mask = self.mask)
        cmd = 'casapy --nologger --log2term -c %s' % template_imager
        exec_cmd(cmd)

        if clean: os.system('rm -rf %s.mask %s.flux %s.psf' % (image, image, image))

    def get_results(self):
        """
        Return image and model name
        """
        # return globbing so to catch any possible tt# for nterm>1
        return glob.glob('%s.image*' % self.image), glob.glob('%s.model*' % self.image)
