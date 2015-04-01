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


class MakeImage(Action):
    """
    Action to make an image

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
    def __init__(self, op_parset, vis_datamap, p, mask_datamap=None,
    	prefix=None, direction=None, band=None, clean=True, index=None,
    	name='MakeImage'):
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
        mask_datamap : Datamap, optional
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
        super(MakeImage, self).__init__(op_parset, name, prefix=prefix,
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
            try:
                os.makedirs(self.image_dir)
            except OSError:
                pass
        self.working_dir = self.image_dir

        # Define names for output images
        imagebasenames = make_image_basename(self.vis_datamap,
            direction=self.direction, band=self.band, prefix=self.prefix)
        self.imagebasenames = [self.image_dir+bn for bn in imagebasenames]

        # Define imaging parameters
        self.set_imaging_parameters()

        # Set a timeout for casapy runs
        if self.op_parset['imager'].lower() == 'casapy':
            self.timeout = 7200

        # Set up all required files
        self.setup()


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import read_mapfile

        # Make first imaging run data maps:
        #     - input is list of MS files
        #     - output is list of image names
        self.p['vis_datamap'] = self.vis_datamap
        self.p['mask_datamap'] = self.mask_datamap

        vis_files, vis_hosts = read_mapfile(self.vis_datamap)
        if self.p['image_final']:
            imagebasenames = self.imagebasenames + '_final'
        else:
            imagebasenames = self.imagebasenames
        self.p['output_datamap'] = self.write_mapfile(imagebasenames,
            prefix=self.prefix+'-imager_output', index=self.index,
            direction=self.direction, band=self.band, host_list=vis_hosts)


    def set_imaging_parameters(self):
        """
        Set various imaging parameters used for both casapy runs
        """
        self.p['cycfactor'] = 3.0
        if self.p['niter'] > 1500:
            self.p['cycfactor'] = 4.0
        self.p['scales'] = []
        if self.p['mscale']:
            self.p['nscales'] = 6
            self.p['scales'] = [0, 3, 7, 25, 60, 150]
            self.p['wsclean_multiscale'] = '' #'-multiscale, ' # need comma for parset
        else:
            self.p['nscales'] = 1
            self.p['scales'] = [0]
            self.p['wsclean_multiscale'] = ''
        self.p['timer'] = ''
        self.p['nfacets'] = 1
        self.imsize = self.p['imsize']
        if self.direction is not None:
            if hasattr(self.direction, 'imsize'):
                self.imsize = self.direction.imsize
        self.p['imsize'] = getOptimumSize(int(self.imsize))
        if 'wplanes' not in self.p:
            # calculate wplanes assuming 1.5" cellsize
            self.p['wplanes'] = 1
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
        if self.op_parset['use_chgcentre']:
            # Let WSClean decide the best value for itself
            self.p['wsclean_wplanes'] = ''
        else:
            self.p['wsclean_wplanes'] = '-nwlayers, {0},'.format(self.p['wplanes'])


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        from factor.lib.datamap_lib import read_mapfile
        from factor.lib.action_lib import get_val_from_str

        if 'ncpu' not in self.p:
            self.p['ncpu'] = self.max_cpu
        if 'n_per_node' not in self.p:
            self.p['n_per_node'] = self.max_cpu

        # Set up mask
        # TODO: use mask datamap
        if self.op_parset['imager'].lower() == 'wsclean':
            if self.mask_datamap is not None:
                mask_file, _ = read_mapfile(self.mask_datamap)
                self.p['mask'] = '-fitsmask, {0}, '.format(mask_file[0])
            else:
                self.p['mask'] = ''
            # TODO: deal with region mask
        elif self.op_parset['imager'].lower() == 'awimager':
            if self.mask_datamap is not None:
                mask_file, _ = read_mapfile(self.mask_datamap)
                self.p['mask'] = '{0}'.format(mask_file[0])
            else:
                self.p['mask'] = ''
        else:
            if self.mask_datamap is None and self.direction is None:
                self.p['mask'] = ['']
            else:
                self.p['mask'] = []
            if self.mask_datamap is not None:
                mask_file, _ = read_mapfile(self.mask_datamap)
                self.p['mask'] += mask_file
            if self.direction is not None and self.op_parset['imager'].lower() == 'casapy':
                # TODO: deal with region mask and AWimager
                self.p['mask'] += [self.direction.reg]
        self.p['imagerroot'] = self.op_parset['imagerroot']

        if self.op_parset['imager'].lower() == 'awimager':
            template = env.get_template('make_image_awimager.pipeline.parset.tpl')
        elif self.op_parset['imager'].lower() == 'casapy':
            template = env.get_template('make_image_casapy.pipeline.parset.tpl')
        elif self.op_parset['imager'].lower() == 'wsclean':
            template = env.get_template('make_image_wsclean.pipeline.parset.tpl')
            self.p['cell_deg'] = get_val_from_str(self.p['cell'], 'deg')
            self.p['threshold_jy'] = get_val_from_str(self.p['threshold'], 'Jy')
            if self.p['nterms'] > 1 and self.direction is not None:
                # nterms > 1 should only be used with a direction
                self.p['nchannels'] = self.direction.nchannels
            else:
                self.p['nchannels'] = 1
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
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
        Remove unneeded files
        """
        pass


class MakeMask(Action):
    """
    Action to make a clean mask from an image
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None,
        direction=None, band=None, clean=True, index=None):
        super(MakeMask, self).__init__(op_parset, 'MakeMask', prefix=prefix,
            direction=direction, band=band, index=index)

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'make_mask'
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
        self.script_file = self.parsetbasename + 'make_clean_mask.py'

        # Set up all required files
        self.setup()


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import read_mapfile

        # Input is list of image basenames
        # Output is clean mask files
        imagebasenames, hosts = read_mapfile(self.input_datamap)
        if self.op_parset['imager'].lower() == 'wsclean':
            if self.p['nterms'] > 1:
                input_files = [bn+'-MFS-image.fits' for bn in imagebasenames]
            else:
                input_files = [bn+'-image.fits' for bn in imagebasenames]
            output_files = [infile+'.cleanmask.fits' for infile in input_files]
        elif self.op_parset['imager'].lower() == 'awimager':
            if self.p['nterms'] == 1:
                input_files = [bn+'.image.restored' for bn in imagebasenames]
            else:
                input_files = [bn+'.image.tt0.restored' for bn in imagebasenames]
            output_files = [bn+'.cleanmask' for bn in imagebasenames]
        else:
            if self.p['nterms'] == 1:
                input_files = [bn+'.image' for bn in imagebasenames]
            else:
                input_files = [bn+'.image.tt0' for bn in imagebasenames]
            output_files = [bn+'.cleanmask' for bn in imagebasenames]

        self.p['input_datamap'] = self.write_mapfile(input_files,
            prefix=self.prefix+'_input', direction=self.direction,
            index=self.index, band=self.band, host_list=hosts)
        self.p['output_datamap'] = self.write_mapfile(output_files,
            prefix=self.prefix+'_output', direction=self.direction,
            index=self.index, band=self.band, host_list=hosts)


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        from factor.lib.datamap_lib import read_mapfile

        if 'ncpu' not in self.p:
            self.p['ncpu'] = self.max_cpu
        if 'n_per_node' not in self.p:
            self.p['n_per_node'] = self.max_cpu

        if self.op_parset['imager'].lower() == 'wsclean' and self.p['nterms'] > 1:
            # Get beam for WSClean images: since MFS image does not (yet) have
            # the beam in its header, we look in the first channel image instead
            from astropy.io import fits

            image_basenames, _ = read_mapfile(self.input_datamap)
            for bn in image_basenames:
                image_file = bn + '-0000-image.fits'
                fits_file = fits.open(image_file, mode="readonly",
                    ignore_missing_end=True)
                hdr = fits_file[0].header
                self.p['beam'] = '({0}, {1}, {2})'.format(hdr['BMAJ'],
                    hdr['BMIN'], hdr['BPA'])
        else:
            self.p['beam'] = 'None'

        self.p['scriptname'] = os.path.abspath(self.script_file)
        template = env.get_template('make_clean_mask.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        if self.op_parset['imager'].lower() == 'wsclean':
            self.p['format'] = 'fits'
        else:
            self.p['format'] = 'casa'

        template = env.get_template('make_clean_mask.tpl')
        tmp = template.render(self.p)
        with open(self.script_file, 'w') as f:
            f.write(tmp)


    def get_results(self):
        """
        Return mask file names.

        If the action has an associated direction, modify the mask to exclude
        regions outside of the direction facet
        """
        from factor.directions import Polygon
        from factor.lib.datamap_lib import read_mapfile
        import pyrap.images as pim
        import numpy as np

        if self.direction is not None:
            output_files = []
            RAverts = self.direction.vertices[0]
            Decverts = self.direction.vertices[1]

            maskfiles, hosts = read_mapfile(self.p['output_datamap'])
            for maskfile in maskfiles:
                parts = maskfile.split('.cleanmask')
                outfile = parts[0] + '.facet_cleanmask' + parts[1]

                mask_im = pim.image(maskfile)
                img_type = mask_im.imagetype()
                if img_type == 'FITSImage':
                    mask_im.saveas(outfile+'tmp')
                    mask_im = pim.image(outfile+'tmp')

                xvert = []
                yvert = []
                for RAvert, Decvert in zip(RAverts, Decverts):
                    pixels = mask_im.topixel([0, 1, Decvert*np.pi/180.0,
                        RAvert*np.pi/180.0])
                    xvert.append(pixels[2]) # x -> Dec
                    yvert.append(pixels[3]) # y -> RA
                poly = Polygon(xvert, yvert)

                # Find masked regions
                data = mask_im.getdata()
                masked_ind = np.where(data[0, 0])

                # Find distance to nearest poly edge and unmask those that
                # are outside the facet
                dist = poly.is_inside(masked_ind[0], masked_ind[1])
                outside_ind = np.where(dist < 0.0)
                if len(outside_ind[0]) > 0:
                    data[0, 0, masked_ind[0][outside_ind], [masked_ind[1][outside_ind]]] = 0

                    # Save changes
                    output_files.append(outfile)
                    mask_im.putdata(data)
                    if img_type == 'FITSImage':
                        mask_im.tofits(outfile, overwrite=True)
                    else:
                        mask_im.saveas(outfile, overwrite=True)

                    # Copy log file that holds clipped rms
                    os.system('cp {0} {1}'.format(maskfile+'.log', outfile+'.log'))
                else:
                    output_files.append(maskfile)

            self.p['output_datamap'] = self.write_mapfile(output_files,
                prefix=self.prefix+'_output', direction=self.direction,
                index=self.index, band=self.band, host_list=hosts)

        return self.p['output_datamap']


class MakeImageIterate(Action):
    """
    Compound action to make iteratively an image using a clean mask
    """
    def __init__(self, op_parset, vis_datamap, p, prefix=None,
        direction=None, band=None, clean=True, index=None):
        super(MakeImageIterate, self).__init__(op_parset, 'MakeImageIterate',
        	prefix=prefix, direction=direction, band=band, index=index,
        	set_up=False)

        self.vis_datamap = vis_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'make_image_iter'
        self.clean = clean


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        pass


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        pass


    def run(self):
        """
        Runs the compound action
        """
        from factor.lib.datamap_lib import read_mapfile

        vis_datamap = self.vis_datamap
        mask_datamap = None
        threshold_5rms = self.p['threshold']

        if self.op_parset['imager'].lower() == 'wsclean':
            # Can't resume with WSclean, so use all iterations on each
            # pass
            self.p['niter'] *= self.p['ncycles']
            self.p['ncycles'] = 2

        for i in range(self.p['ncycles']):
            if self.p['use_rms'] and i == self.p['ncycles'] - 1:
                self.p['threshold'] = threshold_5rms
                self.p['niter'] = 1000000
                if self.op_parset['imager'].lower() == 'casapy':
                    imager.timeout = 7200 # two hours for final pass

            imager = MakeImage(self.op_parset, vis_datamap, self.p,
            	mask_datamap=mask_datamap, direction=self.direction,
            	prefix=self.prefix, band=self.band, index=i)
            if self.op_parset['imager'].lower() == 'casapy':
                imager.timeout = 1800 # half an hour for iterations
            image_basename_mapfile = imager.run()

            if i == self.p['ncycles'] - 1:
                break

            if i > 0 and self.p['iterate_threshold']:
                # Only iterate the threshold for the first pass
                self.p['iterate_threshold'] = False
            masker = MakeMask(self.op_parset, image_basename_mapfile, self.p,
                prefix=self.prefix, direction=self.direction, band=self.band,
                index=i)
            mask_datamap = masker.run()

            mask_file, _ = read_mapfile(mask_datamap)
            log_file = mask_file[0] + '.log'
            with open(log_file, 'r') as f:
                lines = f.readlines()
                threshold_5rms = lines[0].split(': ')[-1] + 'Jy'

        if self.p['image_final']:
            if self.p['use_rms']:
                self.p['threshold'] = threshold_5rms
                self.p['niter'] = 1000000
            imager = MakeImage(self.op_parset, vis_datamap, self.p,
            	mask_datamap=mask_datamap, direction=self.direction,
            	prefix=self.prefix+'_final', band=self.band)
            if self.op_parset['imager'].lower() == 'casapy':
                imager.timeout = 7200 # two hours for final pass
            image_basename_mapfile = imager.run()

        self.output_datamap = image_basename_mapfile

        return self.get_results()



