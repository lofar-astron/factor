"""
Module that holds all image-related actions

Classes
-------
Casapy : Action
    Makes an image with CASA
MakeImage : Action
    Make an image using a clean mask
ExpandMask : Action
    Expand a clean mask to match a larger image

"""

from factor.lib.action import Action
from factor.lib.action_lib import make_basename, make_image_basename
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class Casapy(Action):
    """
    Action to make an image with casapy clean()

    Input data maps
    ---------------
    vis_map : Datamap
        Map of MS files
    mask_map : Datamap, optional
        Map of masks

    Output data maps
    ----------------
    image_map : Datamap
        Map of image basenames

    """
    def __init__(self, op_parset, vis_datamap, p, mask_datamap=None, prefix=None,
        direction=None, localdir=None, clean=True, index=None, image_twice=True,
        name='Casapy'):
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
            Input data map of mask images
        prefix : str, optional
            Prefix to use for image names
        direction : Direction object, optional
            Direction for this image
        localdir : str, optional
            Path to local node directory to use during imaging
        clean : bool, optional
            Remove unneeded files?
        index : int, optional
            Index of action
        image_twice : bool, optional
            If True, image two times, making a clean mask in between. If
            False, image only once; in this case, no clean mask is made.

        """
        super(Casapy, self).__init__(op_parset, name, prefix=prefix,
            direction=direction, index=index)

        # Store input parameters
        self.vis_datamap = vis_datamap
        self.mask_datamap = mask_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'make_image'
        self.localdir = localdir
        self.clean = clean
        self.image_twice = image_twice
        factor_working_dir = op_parset['dir_working']
        if self.direction is None:
            self.image_dir = '{0}/images/{1}/'.format(factor_working_dir,
                self.op_name)
        else:
            self.image_dir = '{0}/images/{1}/{2}/'.format(factor_working_dir,
                self.op_name, self.direction)
        if not os.path.exists(self.image_dir):
            os.makedirs(self.image_dir)
        self.working_dir = self.image_dir

        # Define mask script name
        self.mask_script_file = self.parsetbasename + 'make_clean_mask.py'

        # Define names for output images
        imagebasenames1 = make_image_basename(self.vis_datamap,
            direction=self.direction, prefix=self.prefix+'1')
        self.imagebasenames1 = [self.image_dir+bn for bn in imagebasenames1]
        imagebasenames2 = make_image_basename(self.vis_datamap,
            direction=self.direction, prefix=self.prefix+'2')
        self.imagebasenames2 = [self.image_dir+bn for bn in imagebasenames2]

        # Define imaging parameters
        self.set_imaging_parameters()


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import write_mapfile

        # Make first imaging run data maps:
        #     - input is list of MS files
        #     - output is list of image names
        self.p['vis_datamap_image1'] = self.vis_datamap
        if self.mask_datamap is not None:
            self.p['mask_datamap_image1'] = self.mask_datamap
        self.p['output_datamap_image1'] = write_mapfile(self.imagebasenames1,
            self.op_name, self.name, prefix=self.prefix+'-imager1_output',
            direction=self.direction)

        # Make masking run data maps:
        #     - input is list of images
        #     - output is none
        imnames = []
        masknames = []
        for bn in self.imagebasenames1:
            if self.p['nterms'] > 1:
                imnames.append(bn+'.image.tt0')
            else:
                imnames.append(bn+'.image')
            masknames.append(bn+'.cleanmask')
        self.p['input_datamap_mask'] = write_mapfile(imnames, self.op_name,
            self.name, prefix=self.prefix+'-masker_input', direction=self.direction,
            working_dir=self.op_parset['dir_working'])
        self.p['output_datamap_mask'] = write_mapfile(masknames, self.op_name,
            self.name, prefix=self.prefix+'-masker_output', direction=self.direction,
            working_dir=self.op_parset['dir_working'])

        # Make second imaging run data maps
        #     - input is list of MS files (same as imager 1 inputs)
        #     - output is list of image basenames
        self.p['output_datamap_image2'] = write_mapfile(self.imagebasenames2,
            self.op_name, self.name, prefix=self.prefix+'-imager2_output',
            direction=self.direction)


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
        self.p['imsize'] = self.getOptimumSize(int(self.p['imsize']))
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
        self.p['maskscriptname'] = os.path.abspath(self.mask_script_file)
        if self.image_twice:
            if self.mask_datamap is not None:
                template = env.get_template('make_image_masked.pipeline.parset.tpl')
            else:
                template = env.get_template('make_image.pipeline.parset.tpl')
        else:
            template = env.get_template('make_image_single.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        template_mask = env.get_template('make_clean_mask.tpl')
        tmp = template_mask.render(self.p)
        with open(self.mask_script_file, 'w') as f:
            f.write(tmp)


    def make_pipeline_config_parset(self):
        """
        Writes the pipeline configuration parset
        """
        template = env.get_template('pipeline.cfg.tpl')
        tmp = template.render(self.op_parset)
        with open(self.pipeline_config_file, 'w') as f:
            f.write(tmp)


    def get_results(self):
        """
        Makes and returns map files for image basenames

        Returns
        -------
        images_mapfile : Data map filename
            Filename to map file with imagebasenames

        """
        return self.p['output_datamap_image2']


    def clean(self):
        """
        Remove unneeded files:

        E.g., images/InitSubtract/init_highres/init_highrespipe directory?

        """
        pass


    def getOptimumSize(self, size):
        """
        Gets the nearest optimum image size

        Taken from the casa source code (cleanhelper.py)

        Parameters
        ----------
        size : int
            Target image size in pixels

        Returns
        -------
        optimum_size : int
            Optimum image size nearest to target size

        """
        import numpy

        def prime_factors(n, douniq=True):
            """ Return the prime factors of the given number. """
            factors = []
            lastresult = n
            sqlast=int(numpy.sqrt(n))+1
            if n == 1:
                return [1]
            c=2
            while 1:
                 if (lastresult == 1) or (c > sqlast):
                     break
                 sqlast=int(numpy.sqrt(lastresult))+1
                 while 1:
                     if(c > sqlast):
                         c=lastresult
                         break
                     if lastresult % c == 0:
                         break
                     c += 1

                 factors.append(c)
                 lastresult /= c

            if (factors==[]): factors=[n]
            return  numpy.unique(factors).tolist() if douniq else factors

        n = int(size)
        if (n%2 != 0):
            n+=1
        fac=prime_factors(n, False)
        for k in range(len(fac)):
            if (fac[k] > 7):
                val=fac[k]
                while (numpy.max(prime_factors(val)) > 7):
                    val +=1
                fac[k]=val
        newlarge=numpy.product(fac)
        for k in range(n, newlarge, 2):
            if ((numpy.max(prime_factors(k)) < 8)):
                return k
        return newlarge


class MakeImage(Casapy):
    """
    Action to make an image using a clean mask
    """
    def __init__(self, op_parset, vis_datamap, p, mask_datamap=None, prefix=None,
        direction=None, localdir=None, clean=True, index=None, image_twice=True):
        super(MakeImage, self).__init__(op_parset, vis_datamap, p, prefix=prefix,
            direction=direction, localdir=localdir, clean=clean, index=index,
            image_twice=image_twice, name='MakeImage')
