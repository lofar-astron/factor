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
        if 'timeout' in self.p:
            self.timeout = self.p['timeout']
        else:
            self.timeout = 7200

        # Define names for output images
        imagebasenames = make_image_basename(self.vis_datamap,
            direction=self.direction, band=self.band, prefix=self.prefix)
        self.imagebasenames = [self.image_dir+bn for bn in imagebasenames]
        self.script_file = self.parsetbasename + '_image.py'

        # Define imaging parameters
        self.set_imaging_parameters()

        # Set a timeout for casapy runs
        if self.op_parset['imager'].lower() == 'casapy':
            self.completed_files = []
            self.timeout = 60 # check every 1 minute
            # Define completed file (assumes single image only)
            completed_file = self.imagebasenames[0] + '.done'
            if self.index is not None:
                completed_file += '{0}'.format(self.index)
            self.completed_files.append(completed_file)

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

        # For casapy runs, make datamap for completion files
        if self.op_parset['imager'].lower() == 'casapy':
            self.p['completed_datamap'] = self.write_mapfile(self.completed_files,
                prefix=self.prefix+'-imager_completed', index=self.index,
                direction=self.direction, band=self.band, host_list=vis_hosts)


    def set_imaging_parameters(self):
        """
        Set various imaging parameters
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
            self.p['wsclean_wplanes'] = '' #'-nwlayers, {0},'.format(self.p['wplanes'])


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        from factor.lib.datamap_lib import read_mapfile
        from factor.lib.action_lib import get_val_from_str
        import numpy as np

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
            self.p['scriptname'] = os.path.abspath(self.script_file)
            template = env.get_template('make_image_casapy_comp.pipeline.parset.tpl')
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

        if self.op_parset['imager'].lower() == 'casapy':
            # For casapy only, make a clean script
            template = env.get_template('casapy_clean.tpl')
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
        if 'mask_border' in self.p:
            self.mask_border = self.p['mask_border']
        else:
            self.mask_border = False
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
        import numpy as np

        # Input is list of image basenames
        # Output is clean mask files
        imagebasenames, hosts = read_mapfile(self.input_datamap)
        if self.op_parset['imager'].lower() == 'wsclean':
            if self.p['nterms'] > 1 and self.direction is not None:
                if self.direction.nchannels > 1:
                    input_files = [bn+'-MFS-image.fits' for bn in imagebasenames]
                else:
                    input_files = [bn+'-image.fits' for bn in imagebasenames]
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

        if self.op_parset['imager'].lower() == 'wsclean' and self.direction is not None:
            if self.direction.nchannels > 1:
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

        if self.mask_border:
            output_files = []

            maskfiles, hosts = read_mapfile(self.p['output_datamap'])
            for maskfile in maskfiles:
                parts = maskfile.split('.cleanmask')
                outfile = parts[0] + '.border_cleanmask' + parts[1]

                mask_im = pim.image(maskfile)
                img_type = mask_im.imagetype()
                if img_type == 'FITSImage':
                    mask_im.saveas(outfile+'tmp')
                    mask_im = pim.image(outfile+'tmp')

                # Find masked regions
                data = mask_im.getdata()
                masked_ind = np.where(data[0, 0])

                # Mask pixels along the border
                sh = np.shape(data)
                edge = 25
                data[0, 0, 0:sh[2], 0:edge] = 0
                data[0, 0, 0:edge, 0:sh[3]] = 0
                data[0, 0, 0:sh[2], sh[3]-edge:sh[3]] = 0
                data[0, 0, sh[2]-edge:sh[2], 0:sh[3]] = 0

                # Save changes
                output_files.append(outfile)
                mask_im.putdata(data)
                if img_type == 'FITSImage':
                    mask_im.tofits(outfile, overwrite=True)
                else:
                    mask_im.saveas(outfile, overwrite=True)

                # Copy log file that holds clipped rms
                os.system('cp {0} {1}'.format(maskfile+'.log', outfile+'.log'))

            self.p['output_datamap'] = self.write_mapfile(output_files,
                prefix=self.prefix+'_output', direction=self.direction,
                index=self.index, band=self.band, host_list=hosts)

        if self.direction is not None:
            output_files = []
            RAverts = self.direction.vertices[0]
            Decverts = self.direction.vertices[1]

            maskfiles, hosts = read_mapfile(self.p['output_datamap'])
            for maskfile in maskfiles:
                parts = maskfile.split('cleanmask')
                outfile = parts[0] + 'facet_cleanmask' + parts[1]

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
                # are outside the facet (dist < 0)
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


def image_with_mask(op, imager_parset, prefix, input_mapfiles, directions=None,
    bands=None):
    """
    Helper function to run imaging with masking
    """
    from factor.actions.images import MakeImage, MakeMask
    from factor.lib.datamap_lib import read_mapfile

    mask_datamaps = []
    if directions is not None:
        for direction in directions:
            if direction.region_selfcal is not None and \
                op.parset['imager'].lower() == 'casapy':
                # Use clean region file for first image
                mask_datamaps.append(op.write_mapfile([direction.region_selfcal],
                    prefix=prefix+'_input', direction=direction,
                    index=index))
            else:
                mask_datamaps.append(None)
    else:
        mask_datamaps = [None] * len(input_mapfiles)
    threshold_5rms = [imager_parset['threshold']] * len(input_mapfiles)

    if directions is None:
        directions = [None] * len(input_mapfiles)
    if bands is None:
        bands = [None] * len(input_mapfiles)
    imager_parsets = [imager_parset.copy()] * len(input_mapfiles)

    for i in range(imager_parset['ncycles']):
        if imager_parset['use_rms'] and i == imager_parset['ncycles'] - 1:
            for j, thresh in enumerate(threshold_5rms):
                imager_parsets[j]['threshold'] = thresh
                imager_parsets[j]['niter'] = 1000000

        actions = [MakeImage(op.parset, dm, p, mask_datamap=mm,
            prefix=prefix, direction=d, index=i) for d, dm, mm, p in
            zip(directions, input_mapfiles, mask_datamaps, imager_parsets)]
        image_basename_mapfiles = op.s.run(actions)

        if i == imager_parset['ncycles'] - 1:
            break

        if i > 0 and imager_parset['iterate_threshold']:
            # Only iterate the threshold for the first pass
            imager_parset['iterate_threshold'] = False
        actions = [MakeMask(op.parset, dm, imager_parset,
            prefix=prefix, direction=d, index=i) for d, dm in
            zip(directions, image_basename_mapfiles)]
        mask_mapfiles = op.s.run(actions)

        threshold_5rms = []
        for mm in mask_mapfiles:
            mask_file, _ = read_mapfile(mm)
            log_file = mask_file[0] + '.log'
            with open(log_file, 'r') as f:
                lines = f.readlines()
                threshold_5rms.append(lines[0].split(': ')[-1] + 'Jy')

    return image_basename_mapfiles


