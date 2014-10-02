"""
Action: imager
Make an image
"""

import os
from factor.lib.action import action
from factor.lib.action_lib import makeimagename

from jinja2 import Environment, FileSystemLoader
import os


DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class imager(action):
    """
    Implement the imager action
    """

    def __init__(self, op_name, ms, prefix='', niter=None, imsize=None, cell=None, uvrange=None, mask=None, clean=True):
        super(imager, self).__init__(op_name, name = 'imager')
        self.d = {"ms": ms,
                  "prefix": prefix,
                  "niter": niter,
                  "imsize": imsize,
                  "cell": cell,
                  "uvrange": uvrange,
                  "mask": mask,
                  "clean": clean, # Do we need this here?
                  "image": makeimagename(ms, prefix),
                  }
        self.clean = clean
    
    def _get_command(self):
        template_imager = env.get_template('imager.tpl')
        cmd = 'casapy --nologger --log2term -c %s' % template_imager.render(self.d)
    
    def run(self):
        exec_cmd(self.cmd)

        if clean: 
            os.system('rm -rf %s.mask %s.flux %s.psf' % (image, image, image))

    def get_results(self):
        """
        Return image and model name
        """
        # return globbing so to catch any possible tt# for nterm>1
        return glob.glob('%s.image*' % self.image), glob.glob('%s.model*' % self.image)
