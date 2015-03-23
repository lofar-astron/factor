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
        from factor.lib.operation_lib import make_chunks, merge_chunks

        bands = self.bands

        # Check operation state
        if os.path.exists(self.statebasename+'.done'):
            merged_skymodels_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/InitSubtract/MergeSkymodels/merge_output.datamap')
            skymodels, _ = read_mapfile(merged_skymodels_mapfile)
            for band, skymodel in zip(bands, skymodels):
                band.skymodel_dirindep = skymodel
            merged_data_mapfiles = os.path.join(self.parset['dir_working'],
                'datamaps/InitSubtract/Split/lowres_merged_vis.datamap')
            for band, merged_data_mapfile in zip(bands, merged_data_mapfiles):
                f, _ = read_mapfile(merged_data_mapfile)
                band.file = f[0]
            return

        # Make initial data maps for the input datasets and their dir-indep
        # instrument parmdbs.
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
        actions = [Split(self.parset, dm, p['split'],
            prefix='highres', band=band) for dm, band in zip(vis_mapfiles, bands)]
        split_files_mapfiles= self.s.run(actions)

        self.log.info('High-res imaging...')
        actions = [MakeImageIterate(self.parset, dm, p['imagerh'],
            prefix='highres', band=band) for dm, band in
            zip(split_files_mapfiles, bands)]
        highres_image_basenames_mapfiles = self.s_imager.run(actions)
        basenames = []
        hosts = []
        for bm in highres_image_basenames_mapfiles:
            file_list, host_list = read_mapfile(bm)
            basenames += file_list
            hosts += host_list
        highres_image_basenames_mapfile = self.write_mapfile(basenames,
        	prefix='highres_basenames', host_list=hosts)

        self.log.info('Making high-res sky model...')
        action = MakeSkymodelFromModelImage(self.parset,
            highres_image_basenames_mapfile, p['modelh'], prefix='highres')
        highres_skymodels_mapfile = self.s.run(action)

        if self.parset['use_ftw']:
            self.log.debug('FFTing high-res model image...')
            p['modelh']['imsize'] = p['imagerh']['imsize']
            action = FFT(self.parset, input_data_mapfile,
                highres_image_basenames_mapfile, p['modelh'], prefix='highres')
            self.s.run(action)
            highres_skymodels_mapfile = None

        self.log.debug('Dividing dataset into chunks...')
        chunks_list = []
        for band in bands:
            total_time = (band.endtime - band.starttime) / 3600.0 # hours
            ncpus = self.parset['cluster_specific']['ncpu'] * len(self.parset['cluster_specific']['node_list'])
            chunk_time = np.ceil(total_time/ncpus) + 0.01
            chunks_list.append(make_chunks(band.file, chunk_time,
            	self.parset, 'initsub_chunk', clobber=True))
        chunk_data_mapfiles = []
        chunk_parmdb_mapfiles = []
        chunk_model_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_data_mapfiles.append(self.write_mapfile([chunk.file for chunk in chunks],
                prefix='chunks_vis', direction=d_list[i], host_list=d_hosts[i]))
            if self.parset['use_ftw']:
                chunk_model_mapfiles.append(None)
            else:
                skymodel, hosts = read_mapfile(highres_skymodels_mapfile)
                chunk_model_mapfiles.append(self.write_mapfile([skymodel[i]]*len(chunks),
                    prefix='chunks_highres_skymodel', host_list=hosts[i]))
            parmdb_file, hosts = read_mapfile(dir_indep_parmdbs_mapfile)
            chunk_parmdb_mapfiles.append(self.write_mapfile([parmdb_file[i]]*len(chunks),
                prefix='chunk_parmdb', host_list=hosts[i]))

        self.log.info('Subtracting high-res sky model...')
        actions = [Subtract(self.parset, dm, p['calibh'],
            mm, pm, prefix='highres') for dm, mm, pm in zip(chunk_data_mapfiles,
            chunk_model_mapfiles, chunk_parmdb_mapfiles)]
        self.s.run(actions)

        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(
            	[chunk.file for chunk in chunks], prefix='highres_merge')],
            	prefix='highres_merged_vis', host_list=hosts[i]))

        self.log.info('Averaging...')
        actions = [Average(self.parset, dm, p['avgl'], prefix='highres',
        	band=band) for dm, band in zip(merged_data_mapfiles, bands)]
        avg_files_mapfiles = self.s.run(actions)

        self.log.info('Low-res imaging...')
        actions = [MakeImageIterate(self.parset, dm, p['imagerl'],
            prefix='lowres', band=band) for dm, band in zip(avg_files_mapfiles, bands)]
        lowres_image_basenames_mapfiles = self.s_imager.run(actions)
        basenames = []
        hosts = []
        for bm in lowres_image_basenames_mapfiles:
            file_list, host_list = read_mapfile(bm)
            basenames += file_list
            hosts += host_list
        lowres_image_basenames_mapfile = self.write_mapfile(basenames,
        	prefix='lowres_basenames', host_list=hosts)

        self.log.info('Making low-res sky model...')
        action = MakeSkymodelFromModelImage(self.parset, lowres_image_basenames_mapfile,
            p['modell'], prefix='lowres')
        lowres_skymodels_mapfile = self.s.run(action)
        chunk_model_mapfiles = []
        skymodel, hosts = read_mapfile(lowres_skymodels_mapfile)
        for i, chunks in enumerate(chunks_list):
            if self.parset['use_ftw']:
                chunk_model_mapfiles.append(None)
            else:
                chunk_model_mapfiles.append(self.write_mapfile([skymodel[i]]*
                	len(chunks), prefix='chunks_lowres_skymodel',
                	host_list=hosts[i]))

        if self.parset['use_ftw']:
            self.log.debug('FFTing low-res model image...')
            p['modell']['imsize'] = p['imagerl']['imsize']
            action = [FFT(self.parset, dm, bm, p['modell'], prefix='lowres')
            	for dm, bm in zip(chunk_data_mapfiles,
        	chunk_model_mapfiles)]
            self.s.run(action)
            lowhres_skymodels_mapfile = None

        self.log.info('Subtracting low-res sky model...')
        actions = [Subtract(self.parset, dm, p['calibl'], mm, pm,
        	prefix='highres') for dm, mm, pm in zip(chunk_data_mapfiles,
        	chunk_model_mapfiles, chunk_parmdb_mapfiles)]
        self.s.run(actions)

        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(
                [chunk.file for chunk in chunks], prefix='lowres_merge')],
                prefix='lowres_merged_vis', host_list=hosts[i]))
        for band, merged_data_mapfile in zip(bands, merged_data_mapfiles):
            f, _ = read_mapfile(merged_data_mapfile)
            band.file = f[0]

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
