"""
Action: imager
Make an image
Used parameters:
* niter
* imsize
* cell
* uvrange
* nterms
Return:
* image filename (in a vector to deal with nterm>1)
* model filename (in a vector to deal with nterm>1)
"""

import os
from factor.lib.action import action
from factor.lib.action_lib import makeimagebasename

from jinja2 import Environment, FileSystemLoader
import os


DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class imager(action):
    """
    Implement the imager action
    """

    def __init__(self, op_name, ms, p, prefix = None, direction = None, clean=True):
        super(imager, self).__init__(op_name, name = 'imager')
        self.ms = ms
        self.p = p
        self.imagebasename = makeimagebasename(self.ms, prefix, direction)
        self.clean = clean
    
    def get_command(self):
        # TODO: implement the template
        template_imager = env.get_template('imager.tpl')
        self.cmd = 'casapy --nologger --log2term -c %s' % template_imager.render(self.p)
   
    def run(self):
        self.get_command()
        self.exec_cmd(self.cmd)

        if clean: 
            os.system('rm -rf %s.mask %s.flux %s.psf' % (self.imagebasename, self.imagebasename, self.imagebasename))

    def get_results(self):
        """
        Return image and model name in lists
        """
        # return globbing so to catch any possible tt# for nterm>1
        return glob.glob('%s.image*' % self.image), glob.glob('%s.model*' % self.image)
