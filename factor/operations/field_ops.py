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
        import numpy as np

        bands = self.bands

        # Check operation state
        if os.path.exists(self.statebasename+'.done'):
            merged_skymodels_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/InitSubtract/MergeSkymodels/merge_output.datamap')
            skymodels, _ = read_mapfile(merged_skymodels_mapfile)
            for band, skymodel in zip(bands, skymodels):
                band.skymodel_dirindep = skymodel
            merged_data_mapfiles = [os.path.join(self.parset['dir_working'],
                'datamaps/InitSubtract/{0}/lowres_merged_vis_{0}.datamap'.format(band.name))
                for band in bands]
            for band, merged_data_mapfile in zip(bands, merged_data_mapfile):
                f, _ = read_mapfile(merged_data_mapfile)
                band.file = f[0]
            return

        # Make initial data maps for the input datasets and their dir-indep
        # instrument parmdbs.
        input_data_mapfile = self.write_mapfile([band.file for band in bands],
        	prefix='input_data')
        input_data_mapfiles = []
        files, hosts = read_mapfile(input_data_mapfile)
        for f, h, b in zip(files, hosts, bands):
            input_data_mapfiles.append(self.write_mapfile([f],
        	prefix='input_data_vis', host_list=[h], band=b, index=1))
        dir_indep_parmdbs_mapfile = self.write_mapfile([band.dirindparmdb for band
        	in bands], prefix='dir_indep_parmdbs')

        self.log.info('High-res imaging...')
        actions = [MakeImageIterate(self.parset, dm, p['imagerh'],
            prefix='highres', band=band) for dm, band in
            zip(input_data_mapfiles, bands)]
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
            highres_image_basenames_mapfile, p['imagerh'], prefix='highres')
        highres_skymodels_mapfile = self.s.run(action)

        if self.parset['use_ftw']:
            self.log.debug('FFTing high-res model image...')
            action = FFT(self.parset, input_data_mapfile,
                highres_image_basenames_mapfile, p['imagerh'], prefix='highres')
            self.s_imager.run(action)

        self.log.debug('Dividing datasets into chunks...')
        chunks_list = []
        ncpus = max(int(self.parset['cluster_specific']['ncpu'] * len(
            self.parset['cluster_specific']['node_list']) / len(bands)), 1)
        for i, band in enumerate(bands):
            files, _ = read_mapfile(input_data_mapfiles[i])
            total_time = (band.endtime - band.starttime) / 3600.0 # hours
            chunk_time = min(np.ceil(total_time/ncpus), 1.0) # max of 1 hour per chunk
            chunk_block = int(np.ceil(chunk_time * 3600.0 / band.timepersample))
            self.log.debug('Using {0} time slots ({1:.1f} hr) per chunk for {2}'.format(
                chunk_block, chunk_time, band.name))
            chunks_list.append(make_chunks(files[0], chunk_block,
            	self.parset, 'initsub_chunk', clobber=False))
        chunk_data_mapfiles = []
        chunk_parmdb_mapfiles = []
        chunk_model_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_data_mapfiles.append(self.write_mapfile([chunk.file for chunk in chunks],
                prefix='chunks_vis', band=bands[i], host_list=[hosts[i]]))
            if self.parset['use_ftw']:
                chunk_model_mapfiles.append(None)
            else:
                skymodel, hosts = read_mapfile(highres_skymodels_mapfile)
                chunk_model_mapfiles.append(self.write_mapfile([skymodel[i]]*len(chunks),
                    prefix='chunks_highres_skymodel', band=bands[i], host_list=[hosts[i]]))
            parmdb_file, hosts = read_mapfile(dir_indep_parmdbs_mapfile)
            chunk_parmdb_mapfiles.append(self.write_mapfile([parmdb_file[i]]*len(chunks),
                prefix='chunks_parmdb', band=bands[i], host_list=[hosts[i]]))

        self.log.info('Subtracting high-res sky model...')
        actions = [Subtract(self.parset, dm, p['calibh'],
            mm, pm, prefix='highres', band=band) for dm, mm, pm, band in
            zip(chunk_data_mapfiles, chunk_model_mapfiles, chunk_parmdb_mapfiles,
            bands)]
        self.s.run(actions)

        self.log.debug('Merging chunks for imaging...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(
            	[chunk.file for chunk in chunks], prefix='highres_merge')],
            	prefix='highres_merged_vis', band=bands[i], host_list=[hosts[i]]))

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
            p['imagerl'], prefix='lowres')
        lowres_skymodels_mapfile = self.s.run(action)
        skymodel, hosts = read_mapfile(lowres_skymodels_mapfile)
        chunk_model_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            if self.parset['use_ftw']:
                chunk_model_mapfiles.append(None)
            else:
                chunk_model_mapfiles.append(self.write_mapfile([skymodel[i]]*
                	len(chunks), prefix='chunks_lowres_skymodel',
                	host_list=[hosts[i]]))

        if self.parset['use_ftw']:
            self.log.debug('FFTing low-res model image...')
            actions = [FFT(self.parset, dm, bm, p['imagerl'], band=band,
            	prefix='lowres') for dm, bm, band in zip(merged_data_mapfiles,
            	lowres_image_basenames_mapfiles, bands)]
            self.s_imager.run(actions)

        self.log.debug('Dividing datasets into chunks...')
        chunks_list = []
        for i, band in enumerate(bands):
            files, _ = read_mapfile(merged_data_mapfiles[i])
            total_time = (band.endtime - band.starttime) / 3600.0 # hours
            chunk_time = min(np.ceil(total_time/ncpus), 1.0) # max of 1 hour per chunk
            chunk_block = int(np.ceil(chunk_time * 3600.0 / band.timepersample))
            self.log.debug('Using {0} time slots ({1:.1f} hr) per chunk for {2}'.format(
                chunk_block, chunk_time, band.name))
            chunks_list.append(make_chunks(files[0], chunk_block,
            	self.parset, 'initsub_chunk', clobber=True))
        chunk_data_mapfiles = []
        chunk_parmdb_mapfiles = []
        chunk_model_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_data_mapfiles.append(self.write_mapfile([chunk.file for chunk in chunks],
                prefix='chunks_vis', band=bands[i], host_list=[hosts[i]]))
            if self.parset['use_ftw']:
                chunk_model_mapfiles.append(None)
            else:
                skymodel, hosts = read_mapfile(lowres_skymodels_mapfile)
                chunk_model_mapfiles.append(self.write_mapfile([skymodel[i]]*len(chunks),
                    prefix='chunks_highres_skymodel', band=bands[i], host_list=[hosts[i]]))
            parmdb_file, hosts = read_mapfile(dir_indep_parmdbs_mapfile)
            chunk_parmdb_mapfiles.append(self.write_mapfile([parmdb_file[i]]*len(chunks),
                prefix='chunks_parmdb', band=bands[i], host_list=[hosts[i]]))

        self.log.info('Subtracting low-res sky model...')
        actions = [Subtract(self.parset, dm, p['calibl'], mm, pm,
        	prefix='lowres', band=band) for dm, mm, pm, band in
        	zip(chunk_data_mapfiles, chunk_model_mapfiles,
        	chunk_parmdb_mapfiles, bands)]
        self.s.run(actions)

        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(
                [chunk.file for chunk in chunks], prefix='lowres_merge')],
                prefix='lowres_merged_vis', band=bands[i], host_list=[hosts[i]]))
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
