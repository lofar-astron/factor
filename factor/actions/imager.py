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

    def __init__(self, op_name, ms, p, clean=True):
        super(imager, self).__init__(op_name, name = 'imager')
        self.ms = ms
        self.p = p
        self.clean = clean
        self.image = makeimagename(ms, prefix)

    def run(self):
        # TODO: implement the template
        template_imager = make_template(ms = self.ms, image = self.image, p = self.p)
        cmd = 'casapy --nologger --log2term -c %s' % template_imager
        exec_cmd(cmd)

        if clean: os.system('rm -rf %s.mask %s.flux %s.psf' % (image, image, image))

    def get_results(self):
        """
        Return image and model name in lists
        """
        # return globbing so to catch any possible tt# for nterm>1
        return glob.glob('%s.image*' % self.image), glob.glob('%s.model*' % self.image)
