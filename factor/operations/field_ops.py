"""
Module that holds all field (non-facet-specific) operations

Classes
-------
InitSubtract : Operation
    Images each band at high and low resolution to make and subtract sky models
MakeMosaic : Operation
    Makes a mosaic from the facet images

"""
import os
from factor.lib.operation import Operation
from factor.lib.scheduler import Scheduler


class InitSubtract(Operation):
    """
    Operation to create empty datasets
    """
    def __init__(self, parset, bands, reset=False):
        super(InitSubtract, self).__init__(parset, bands, direction=None,
            reset=reset, name='InitSubtract')

        # Set up imager scheduler (runs at most num_nodes imagers in parallel)
        num_nodes = len(self.parset['cluster_specific']['node_list'])
        self.s_imager = Scheduler(max_threads=num_nodes, name=self.name,
            op_parset=self.parset)


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.images import MakeImageIterate
        from factor.actions.models import MakeSkymodelFromModelImage, MergeSkymodels, FFT
        from factor.actions.calibrations import Subtract
        from factor.actions.visibilities import Average, Split
        from factor.lib.datamap_lib import read_mapfile
        from factor.operations.hardcoded_param import init_subtract as p

        bands = self.bands

        # Check operation state
        if os.path.exists(self.statebasename+'.done'):
            merged_skymodels_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/InitSubtract/MergeSkymodels/merge_output.datamap')
            skymodels, _ = read_mapfile(merged_skymodels_mapfile)
            for band, skymodel in zip(bands, skymodels):
                band.skymodel_dirindep = skymodel
            return

        # Make initial data maps for the input datasets and their dir-indep
        # instrument parmdbs. Also make datamaps with one file each for casapy
        # imaging run (otherwise we get conflicts with temp files)
        input_data_mapfile = self.write_mapfile([band.file for band in bands],
        	prefix='input_data')
        vis_mapfiles = []
        files, hosts = read_mapfile(input_data_mapfile)
        for f, h, b in zip(files, hosts, bands):
            vis_mapfiles.append(self.write_mapfile([f],
        	prefix='input_data_vis', host_list=[h], band=b, index=1))
        dir_indep_parmdbs_mapfile = self.write_mapfile([band.dirindparmdb for band
        	in bands], prefix='dir_indep_parmdbs')

        self.log.info('Spliting off corrected data...')
        action = Split(self.parset, input_data_mapfile, p['split'],
            prefix='highres')
        split_files_mapfile = self.s.run(action)
        vis_mapfiles = []
        files, hosts = read_mapfile(split_files_mapfile)
        for f, h, b in zip(files, hosts, bands):
            vis_mapfiles.append(self.write_mapfile([f],
        	prefix='corrected_vis', host_list=[h], band=b))

        self.log.info('High-res imaging...')
        actions = [MakeImageIterate(self.parset, dm, p['imagerh'],
            prefix='highres', band=band) for dm, band in zip(vis_mapfiles, bands)]
        highres_image_basenames_mapfiles = self.s_imager.run(actions)
        basenames = []
        hosts = []
        for bm in highres_image_basenames_mapfiles:
            file_list, host_list = read_mapfile(bm)
            basenames += file_list
            hosts += host_list
        highres_image_basenames_mapfile = self.write_mapfile(basenames,
        	prefix='highres_basenames', host_list=hosts)

        if self.parset['use_ftw']:
            self.log.debug('FFTing high-res model image...')
            p['modelh']['imsize'] = p['imagerh']['imsize']
            action = FFT(self.parset, input_data_mapfile,
                highres_image_basenames_mapfile, p['modelh'], prefix='highres')
            self.s.run(action)
            highres_skymodels_mapfile = None

        self.log.info('Making high-res sky model...')
        action = MakeSkymodelFromModelImage(self.parset,
            highres_image_basenames_mapfile, p['modelh'], prefix='highres')
        highres_skymodels_mapfile = self.s.run(action)

        self.log.info('Subtracting high-res sky model...')
        action = Subtract(self.parset, input_data_mapfile, p['calibh'],
            model_datamap=highres_skymodels_mapfile,
            parmdb_datamap=dir_indep_parmdbs_mapfile, prefix='highres')
        self.s.run(action)

        self.log.info('Averaging...')
        action = Average(self.parset, input_data_mapfile, p['avgl'],
            prefix='highres')
        avg_files_mapfile = self.s.run(action)
        vis_mapfiles = []
        files, hosts = read_mapfile(avg_files_mapfile)
        for f, h, b in zip(files, hosts, bands):
            vis_mapfiles.append(self.write_mapfile([f],
        	prefix='input_data_vis', host_list=[h], band=b, index=2))

        self.log.info('Low-res imaging...')
        actions = [MakeImageIterate(self.parset, dm, p['imagerl'],
            prefix='lowres', band=band) for dm, band in zip(vis_mapfiles, bands)]
        lowres_image_basenames_mapfiles = self.s_imager.run(actions)
        basenames = []
        hosts = []
        for bm in lowres_image_basenames_mapfiles:
            file_list, host_list = read_mapfile(bm)
            basenames += file_list
            hosts += host_list
        lowres_image_basenames_mapfile = self.write_mapfile(basenames,
        	prefix='lowres_basenames', host_list=hosts)

        if self.parset['use_ftw']:
            self.log.debug('FFTing low-res model image...')
            p['modell']['imsize'] = p['imagerl']['imsize']
            action = FFT(self.parset, input_data_mapfile,
                lowres_image_basenames_mapfile, p['modell'], prefix='lowres')
            self.s.run(action)

        self.log.info('Making low-res sky model...')
        action = MakeSkymodelFromModelImage(self.parset, lowres_image_basenames_mapfile,
            p['modell'], prefix='lowres')
        lowres_skymodels_mapfile = self.s.run(action)

        self.log.info('Subtracting low-res sky model...')
        action = Subtract(self.parset, input_data_mapfile, p['calibl'],
            model_datamap=lowres_skymodels_mapfile,
            parmdb_datamap=dir_indep_parmdbs_mapfile, prefix='lowres')
        self.s.run(action)

        self.log.info('Merging low- and high-res sky models...')
        action = MergeSkymodels(self.parset, lowres_skymodels_mapfile,
            highres_skymodels_mapfile, p['merge'], prefix='merge')
        merged_skymodels_mapfile = self.s.run(action)
        skymodels, _ = read_mapfile(merged_skymodels_mapfile)
        for band, skymodel in zip(bands, skymodels):
            band.skymodel_dirindep = skymodel


class MakeMosaic(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(MakeMosaic, self).__init__(parset, bands, direction=direction,
            reset=reset, name='MakeMosaic')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        pass
