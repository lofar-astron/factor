"""
Module that holds all image-related objects and functions

Actions defined in this module:

- make_image: makes an image using clean masks

"""

from factor.lib.action import Action
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class make_image(Action):
    """
    Makes restored and model images using clean masks

    Returns data map of images and models
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        localdir=None, clean=True, index=None):
        """
        Create make_image action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map of MS file(s) to image
        p : parset dict
            Input parset dict defining imaging and pipeline parameters
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
        """
        from factor.lib.action_lib import make_parset_basename, make_pipeline_dirname

        super(make_image, self).__init__(op_parset, 'make_image')

        self.input_datamap = input_datamap
        self.p = p.copy()
        if prefix is None:
            prefix = 'make_image'
        self.prefix = prefix
        self.direction = direction
        self.localdir = localdir
        self.clean = clean
        self.index = index

        # Set up parset and script names
        self.parsetbasename = self.parset_dir + make_parset_basename(prefix,
            direction, index)
        self.pipeline_parset_file = self.parsetbasename + 'pipe.parset'
        self.pipeline_config_file = self.parsetbasename + 'pipe.cfg'
        self.maskscriptname = self.parsetbasename + 'makecleanmask.py'
        self.pipeline_run_dir += make_pipeline_dirname(prefix, direction, index)
        if not os.path.exists(self.pipeline_run_dir):
            os.makedirs(self.pipeline_run_dir)

        # Set up image names for output data map
        self.imagebasenames1 = []
        self.imagebasenames2 = []
        imagebasenames1 = self.make_image_basename(prefix+'1')
        for bn in imagebasenames1:
            self.imagebasenames1.append(self.image_dir+bn)
        imagebasenames2 = self.make_image_basename(prefix+'2')
        for bn in imagebasenames1:
            self.imagebasenames2.append(self.image_dir+bn)

        # Make data maps
        self.make_datamaps()

        # Make pipeline parsets
        self.set_imaging_parameters()
        self.set_pipeline_parameters()
        self.make_pipeline_control_parset()
        self.make_pipeline_config_parset()

        # Run the pipeline
        self.run()


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import write_mapfile

        # Make first imaging data maps:
        #     - input is list of MS files
        #     - output is list of image names
        self.p['input_datamap_image1'] = self.input_datamap
        self.p['output_datamap_image1'] = write_mapfile(self.imagebasenames1,
            self.op_name, self.name, prefix=self.prefix+'-imager1_output',
            direction=self.direction)

        # Make masking data maps:
        #     - input is list of model images
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
            self.name, prefix=self.prefix+'-masker_input', direction=self.direction)
        self.p['output_datamap_mask'] = write_mapfile(masknames, self.op_name,
            self.name, prefix=self.prefix+'-masker_output', direction=self.direction)

        # Make second imaging data maps
        #     - input is list of MS files (same as imager 1 inputs)
        #     - output is list of image names
        self.p['input_datamap_image2'] = self.input_datamap
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


    def set_pipeline_parameters(self):
        """
        Sets various pipeline parameters specific to this action
        """
        self.op_parset['runtime_dir'] = os.path.join(os.path.abspath('.'), self.pipeline_run_dir)
        self.op_parset['working_dir'] = os.path.join(os.path.abspath('.'), self.image_dir)
        self.op_parset['ncpu'] = 1


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        self.p['maskscriptname'] = os.path.abspath(self.maskscriptname)
        template = env.get_template('imager.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        template_mask = env.get_template('makecleanmask.tpl')
        tmp = template_mask.render(self.p)
        with open(self.maskscriptname, 'w') as f:
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
        Return map files for final images, models, and clean masks

        If nterms > 1, only the tt0 image and model are returned
        """
        from factor.lib.datamap_lib import write_mapfile

        images = []
        models = []
        masks = []
        nterms = self.p['nterms']
        imagebasenames = self.imagebasenames2
        for ib in imagebasenames:
            if nterms > 1:
                images.append(ib+'.image.tt0')
                models.append(ib+'.model.tt0')
            else:
                images.append(ib+'.image')
                models.append(ib+'.model')
            masks.append(ib+'.cleanmask')
        images_mapfile = write_mapfile(images, self.op_name, action_name=self.name,
            prefix=self.prefix, direction=self.direction, index=self.index)
        models_mapfile = write_mapfile(models, self.op_name, action_name=self.name,
            prefix=self.prefix, direction=self.direction, index=self.index)
        masks_mapfile = write_mapfile(masks, self.op_name, action_name=self.name,
            prefix=self.prefix, direction=self.direction, index=self.index)

        return images_mapfile, models_mapfile, masks_mapfile


    def make_image_basename(self, prefix=None):
        """
        define a standard name pattern for imaging files
        """
        from factor.lib.datamap_lib import read_mapfile
        import re
        import os

        if prefix is None:
            prefix = 'image'
            logging.warn('Prefix of imagebasename not selected, using "{0}".'.format(prefix))

        msfiles = read_mapfile(self.input_datamap)
        image_basenames = []

        for msfile in msfiles:
            msbase = os.path.basename(msfile)

            if self.direction is not None:
                try:
                    dirtxt = direction.name
                except:
                    dirtxt = direction
                image_basenames.append('%s-%s_%s' % (prefix, re.sub(r'.MS|.ms', '', msbase), dirtxt))
            else:
                image_basenames.append('%s-%s' % (prefix, re.sub(r'.MS|.ms', '', msbase)))

        return image_basenames


    def getOptimumSize(self, size):
        '''
        Gets the nearest optimal image size

        Taken from the casa source code (cleanhelper.py)
        '''
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

        n=size
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

