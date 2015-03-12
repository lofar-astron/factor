"""
Module that holds all image-related actions

Classes
-------
Casapy : Action
    Makes an image with CASA
MakeImage : Action
    Make an image using a clean mask

"""

from factor.lib.action import Action
from factor.lib.action_lib import make_image_basename, getOptimumSize
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class Casapy(Action):
    """
    Action to make an image with casapy clean()

    Input data maps
    ---------------
    vis_datamap : Datamap
        Map of MS files
    mask_datamap : Datamap, optional
        Map of masks

    Output data maps
    ----------------
    image_datamap : Datamap
        Map of image basenames

    """
    def __init__(self, op_parset, vis_datamap, p, mask_datamap=None, prefix=None,
        direction=None, band=None, clean=True, index=None, name='Casapy'):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        vis_datamap : Datamap
            Input data map of MS file(s) to image
        p : parset dict
            Input parset dict defining imaging and pipeline parameters
        maks_datamap : Datamap, optional
            Input data map of mask images or CASA-format regions
        prefix : str, optional
            Prefix to use for image names
        direction : Direction object, optional
            Direction for this image
        band : Band object, optional
            Band for this image
        clean : bool, optional
            Remove unneeded files?
        index : int, optional
            Index of action

        """
        super(Casapy, self).__init__(op_parset, name, prefix=prefix,
            direction=direction, band=band, index=index)

        # Store input parameters
        self.vis_datamap = vis_datamap
        self.mask_datamap = mask_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'make_image'
        self.clean = clean
        self.image_dir += '{0}/{1}/'.format(self.op_name, self.name)
        if self.direction is not None:
            self.image_dir += '{0}/'.format(self.direction.name)
        if self.band is not None:
            self.image_dir += '{0}/'.format(self.band.name)
        if not os.path.exists(self.image_dir):
            os.makedirs(self.image_dir)
        self.working_dir = self.image_dir

        # Define script name
        self.script_file = self.parsetbasename + 'make_image.py'

        # Define names for output images
        imagebasenames = make_image_basename(self.vis_datamap,
            direction=self.direction, band=self.band, prefix=self.prefix)
        self.imagebasenames = [self.image_dir+bn for bn in imagebasenames]

        # Define imaging parameters
        self.set_imaging_parameters()

        # Set up all required files
        self.setup()


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import read_mapfile

        # Make first imaging-run data maps:
        #     - input is list of MS files
        #     - output is list of image names
        self.p['vis_datamap'] = self.vis_datamap
        self.p['mask_datamap'] = self.mask_datamap

        vis_files, vis_hosts = read_mapfile(self.vis_datamap)
        self.p['output_datamap'] = self.write_mapfile(self.imagebasenames,
            prefix=self.prefix+'-imager_output', index=self.index,
            direction=self.direction, host_list=vis_hosts)


    def set_imaging_parameters(self):
        """
        Set various imaging parameters used for both casapy runs
        """
        self.p['cycfactor'] = 3.0
        if self.p['niter'] > 1500:
            self.p['cycfactor'] = 4.0
        self.p['scales'] = []
        if self.p['mscale']:
            self.p['scales'] = [0,3,7,25,60,150]
        self.p['timer'] = ''
        self.p['nfacets'] = 1
        self.p['wplanes'] = 1
        self.imsize = self.p['imsize']
        if self.direction is not None:
            if hasattr(self.direction, 'imsize'):
                self.imsize = self.direction.imsize
        self.p['imsize'] = getOptimumSize(int(self.imsize))
        if self.p['imsize'] > 512:
            self.p['wplanes'] = 64
        if self.p['imsize'] > 799:
            self.p['wplanes'] = 96
        if self.p['imsize'] > 1023:
            self.p['wplanes'] = 128
        if self.p['imsize'] > 1599:
            self.p['wplanes'] = 256
            self.p['nfacets'] = 1
        if self.p['imsize'] > 2047:
            self.p['wplanes'] = 384
            self.p['nfacets'] = 1
        if self.p['imsize'] > 3000:
            self.p['wplanes'] = 448
            self.p['nfacets'] = 1
        if self.p['imsize'] > 4095:
            self.p['wplanes'] = 512
            self.p['nfacets'] = 1


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        from factor.lib.datamap_lib import read_mapfile

        if 'ncpu' not in self.p:
            self.p['ncpu'] = self.max_cpu
        if self.mask_datamap is None:
            self.p['mask'] = ''
        if self.direction is not None:
            self.p['mask'] = self.direction.reg

        self.p['scriptname'] = os.path.abspath(self.script_file)
        template = env.get_template('make_image.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        template = env.get_template('make_image.tpl')
        tmp = template.render(self.p)
        with open(self.script_file, 'w') as f:
            f.write(tmp)


    def get_results(self):
        """
        Makes and returns map files for image basenames

        Returns
        -------
        images_mapfile : Data map filename
            Filename to map file with imagebasenames

        """
        return self.p['output_datamap']


    def clean(self):
        """
        Remove unneeded files:

        E.g., images/InitSubtract/init_highres/init_highrespipe directory?

        """
        pass


class MakeImage(Casapy):
    """
    Action to make an image using a clean mask
    """
    def __init__(self, op_parset, vis_datamap, p, mask_datamap=None, prefix=None,
        direction=None, band=None, clean=True, index=None):
        super(MakeImage, self).__init__(op_parset, vis_datamap, p, prefix=prefix,
            direction=direction, band=band, clean=clean, index=index,
            name='MakeImage')
