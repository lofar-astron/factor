"""
Module that holds all facet operations

Operations are largely defined by logical groups of actions, but
they also break up the processing into parallel or serial steps.

Classes
-------
FacetAdd : Operation
    Adds facet sources to data
FacetSetup : Operation
    Sets up the data for selfcal
FacetSelfcal : Operation
    Runs selfcal cycle on facet calibrator
FacetImage : Operation
    Images the entire facet
FacetSub : Operation
    Subtracts all facet sources from data and images residuals
FacetAddAllFinal : Operation
    Adds all facet sources from final model
FacetImageFinal : Operation
    Images the entire facet

"""
import os
from factor.lib.operation import Operation
from factor.lib.scheduler import Scheduler


class FacetAdd(Operation):
    """
    Operation to add calibrator source to data
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetAdd, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetAdd')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import PhaseShift
        from factor.actions.models import MakeFacetSkymodel
        from factor.actions.calibrations import Add
        from factor.operations.hardcoded_param import facet_add as p
        from factor.lib.datamap_lib import read_mapfile

        bands = self.bands
        d = self.direction

        # Check state
        if self.check_completed(d):
            return

        if os.path.exists(self.statebasename+'.done'):
            shifted_cal_data_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/FacetAdd/PhaseShift/{0}/facet_output_{0}-1.datamap'.
                format(d.name))
            shifted_all_data_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/FacetAdd/PhaseShift/{0}/facet_output_{0}-2.datamap'.
                format(d.name))
            shifted_cal_data_files, hosts = read_mapfile(shifted_cal_data_mapfile)
            shifted_all_data_files, hosts = read_mapfile(shifted_all_data_mapfile)
            d.shifted_cal_data_files = shifted_cal_data_files
            d.shifted_all_data_files = shifted_all_data_files
            return

        # Make initial data maps for the empty datasets, their dir-indep
        # instrument parmdbs, and their dir-indep sky models
        subtracted_all_mapfile = self.write_mapfile([band.file for band in bands],
        	prefix='subtracted_all', direction=d)
        dir_indep_parmdbs_mapfile = self.write_mapfile([band.dirindparmdb for band
        	in bands], prefix='dir_indep_parmdbs', direction=d)
        dir_indep_skymodels_mapfile = self.write_mapfile([band.skymodel_dirindep
        	for band in bands], prefix='dir_indep_skymodels', direction=d)

        # Add sources from the dir-indep sky model for this direction to the
        # visibilities
        self.log.info('Selecting sources for this direction...')
        action = MakeFacetSkymodel(self.parset, dir_indep_skymodels_mapfile,
            {}, d, prefix='cal', cal_only=True, index=1)
        dir_indep_cal_skymodels_mapfile = self.s.run(action)
        action = MakeFacetSkymodel(self.parset, dir_indep_skymodels_mapfile,
            {}, d, prefix='all', cal_only=False, index=2)
        dir_indep_all_skymodels_mapfile = self.s.run(action)

        self.log.info('Adding sources for this direction...')
        self.parset['use_ftw'] = False
        if bands[0].has_sub_data_new:
            p['add']['incol'] += '_NEW'
        action = Add(self.parset, subtracted_all_mapfile, p['add_cal'],
            dir_indep_cal_skymodels_mapfile, dir_indep_parmdbs_mapfile,
            prefix='facet_dirindep', direction=d, index=1)
        self.s.run(action)
        action = Add(self.parset, subtracted_all_mapfile, p['add_all'],
            dir_indep_all_skymodels_mapfile, dir_indep_parmdbs_mapfile,
            prefix='facet_dirindep', direction=d, index=2)
        self.s.run(action)

        # Phase shift to facet center
        self.log.info('Phase shifting...')
        action = PhaseShift(self.parset, subtracted_all_mapfile, p['shift_cal'],
            prefix='facet', direction=d, index=1)
        shifted_cal_data_mapfile = self.s.run(action)
        shifted_cal_data_files, _ = read_mapfile(shifted_cal_data_mapfile)
        action = PhaseShift(self.parset, subtracted_all_mapfile, p['shift_all'],
            prefix='facet', direction=d, index=2)
        shifted_all_data_mapfile = self.s.run(action)
        shifted_all_data_files, _ = read_mapfile(shifted_all_data_mapfile)
        action = PhaseShift(self.parset, subtracted_all_mapfile, p['shift_sub'],
            prefix='facet', direction=d, index=3)
        shifted_sub_data_mapfile = self.s.run(action)
        shifted_sub_data_files, _ = read_mapfile(shifted_sub_data_mapfile)

        d.shifted_cal_data_files = shifted_cal_data_files
        d.shifted_all_data_files = shifted_all_data_files
        d.shifted_sub_data_files = shifted_sub_data_files

        # Save state
        self.set_completed(d)


class FacetSetup(Operation):
    """
    Operation to set up data for selfcal
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetSetup, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetSetup')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import Average, Concatenate
        from factor.actions.calibrations import Apply
        from factor.lib.operation_lib import copy_column
        from factor.operations.hardcoded_param import facet_setup as p
        from factor.lib.datamap_lib import read_mapfile

        bands = self.bands
        d_list = self.direction
        if type(d_list) is not list:
            d_list = [d_list]

        # Check state
        if self.check_completed(d_list):
            return

        # Divide up the nodes among the directions
        node_list = self.parset['cluster_specific']['node_list']
        if len(d_list) >= len(node_list):
            for i in range(len(d_list)-len(node_list)):
                node_list.append(node_list[i])
            d_hosts = [[n] for n in node_list]
        else:
            parts = len(d_list)
            d_hosts = [node_list[i*len(node_list)//parts:
                (i+1)*len(node_list)//parts] for i in range(parts)]

        # Make initial data maps for the phase-shifted calibrator datasets
        # and their dir-indep instrument parmdbs
        shifted_data_mapfiles = []
        dir_indep_parmdbs_mapfiles = []
        for d, h in zip(d_list, d_hosts):
            shifted_data_mapfiles.append(self.write_mapfile(d.shifted_cal_data_files,
                prefix='shifted', direction=d, host_list=h))
            dir_indep_parmdbs_mapfiles.append(self.write_mapfile([band.
            	dirindparmdb for band in bands], prefix='dir_indep_parmdbs',
            	direction=d, host_list=h))

        # apply direction-independent calibration
        self.log.info('Applying direction-independent calibration...')
        actions = [Apply(self.parset, dm, p['apply'],
            pm, prefix='facet_dirindep', direction=d) for d, dm, pm in zip(d_list,
            shifted_data_mapfiles, dir_indep_parmdbs_mapfiles)]
        self.s.run(actions)

        # average to 1 channel per band. Do this twice, once for DATA and once
        # for CORRECTED_DATA
        self.log.info('Averaging...')
        actions = [Average(self.parset, m, p['avg1'], prefix='facet',
            direction=d, index=1) for d, m in zip(d_list, shifted_data_mapfiles)]
        avg1_data_mapfiles = self.s.run(actions)
        actions = [Average(self.parset, m, p['avg2'], prefix='facet',
            direction=d, index=2) for d, m in zip(d_list, shifted_data_mapfiles)]
        avg2_data_mapfiles = self.s.run(actions)

        # concatenate all phase-shifted, averaged bands together. Do this twice,
        # once for each of the two averaged file sets
        self.log.info('Concatenating bands...')
        actions = [Concatenate(self.parset, m, p['concat1'],
            prefix='facet_bands', direction=d, index=1) for d, m in zip(d_list,
            avg1_data_mapfiles)]
        concat_data_mapfiles = self.s.run(actions)
        actions = [Concatenate(self.parset, m, p['concat2'],
            prefix='facet_bands', direction=d, index=2) for d, m in zip(d_list,
            avg2_data_mapfiles)]
        concat_corrdata_mapfiles = self.s.run(actions)

        # Copy over DATA column (was phase-shifted CORRECTED_DATA) from second
        # concat file, so that first concatenated file has both DATA and
        # dir-indep CORRECTED_DATA
        for dm, cdm, d in zip(concat_data_mapfiles, concat_corrdata_mapfiles,
            d_list):
            concat_data_file, _ = read_mapfile(dm)
            concat_corrdata_file, _ = read_mapfile(cdm)
            d.cal_concat_file = concat_data_file[0]
            copy_column(d.cal_concat_file, p['copy']['incol'], p['copy']['outcol'],
                ms_from=concat_corrdata_file[0])

        # Save state
        self.set_completed(d_list)


class FacetSelfcal(Operation):
    """
    Operation to selfcal one or more directions
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetSelfcal, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetSelfcal')

        # Set up scheduler (runs at most num_nodes directions in parallel)
        num_nodes = len(self.parset['cluster_specific']['node_list'])
        self.s = Scheduler(max_threads=num_nodes, name=self.name,
            op_parset=self.parset)


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import Average, Concatenate, PhaseShift
        from factor.actions.calibrations import Apply, Solve
        from factor.actions.images import MakeImageIterate
        from factor.actions.models import FFT
        from factor.actions.solutions import Smooth, ResetPhases
        from factor.lib.operation_lib import copy_column, make_chunks, \
            merge_chunks, merge_parmdbs, merge_chunk_parmdbs
        from factor.operations.hardcoded_param import facet_selfcal as p
        from factor.lib.datamap_lib import read_mapfile

        bands = self.bands
        d_list = self.direction
        if type(d_list) is not list:
            d_list = [d_list]

        # Check state
        if self.check_completed(d_list):
            return

        # Set imager
        self.parset['imager'] = self.parset['imager_selfcal']

        # Divide up the nodes among the directions
        node_list = self.parset['cluster_specific']['node_list']
        if len(d_list) >= len(node_list):
            for i in range(len(d_list)-len(node_list)):
                node_list.append(node_list[i])
            d_hosts = [[n] for n in node_list]
        else:
            parts = len(d_list)
            d_hosts = [node_list[i*len(node_list)//parts:
                (i+1)*len(node_list)//parts] for i in range(parts)]

        # Make initial data maps for the averaged, phase-shifted datasets
        facet_data_mapfiles = []
        facet_unavg_data_mapfiles = []
        for d, h in zip(d_list, d_hosts):
            facet_data_mapfiles.append(self.write_mapfile([d.cal_concat_file],
                prefix='shifted_vis', direction=d, host_list=h))

        # Set image sizes
        for d in d_list:
            cell = float(p['imager0']['cell'].split('arcsec')[0]) # arcsec per pixel
            imsize = d.cal_radius_deg * 1.5 * 3600.0 / cell # pixels
            if imsize < 512:
                imsize = 512
            d.imsize = imsize

        self.log.info('Imaging (facet image #0)...')
        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg0'], prefix='facet',
            direction=d, index=0) for d, m in zip(d_list, facet_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Imaging...')
        actions = [MakeImageIterate(self.parset, dm, p['imager0'],
        	prefix='facet_selfcal0', direction=d) for d, dm in zip(d_list,
        	avg_data_mapfiles)]
        image0_basenames_mapfiles = self.s.run(actions)

        self.log.info('FFTing model image (facet model #0)...')
        actions = [FFT(self.parset, dm, mm, p['imager0'], prefix='fft0',
            direction=d, index=0) for d, dm, mm in zip(d_list,
            facet_data_mapfiles, image0_basenames_mapfiles)]
        self.s.run(actions)

        self.log.debug('Dividing dataset into chunks...')
        chunks_list = []
        for d, m in zip(d_list, facet_data_mapfiles):
            files, _ = read_mapfile(m)
            chunks_list.append(make_chunks(files[0], d.solint_a*5.0,
            	self.parset, 'facet_chunk', direction=d, clobber=True))
        chunk_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_data_mapfiles.append(self.write_mapfile([chunk.file for chunk in chunks],
                prefix='chunks_vis', direction=d_list[i], host_list=d_hosts[i]))

        self.log.info('Solving for phase solutions and applying them (#1)...')
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseonly1'].copy()
            p_d['timestep'] = d.solint_p
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=None,
            prefix='facet_phaseonly', direction=d, index=0)
            for d, dm, pd in zip(d_list, chunk_data_mapfiles, p_list)]
        self.s.run(actions)

        self.log.info('Imaging (facet image #1)...')
        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg1'], prefix='facet',
            direction=d, index=1) for d, m in zip(d_list, chunk_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Merging averged chunks...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            avg_files, _ = read_mapfile(avg_data_mapfiles[i])
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(avg_files,
            prefix=None, clobber=True)], direction=d_list[i],
            host_list=d_hosts[i]))

        self.log.debug('Imaging...')
        actions = [MakeImageIterate(self.parset, m, p['imager1'],
            prefix='facet_selfcal1', direction=d) for d, m in zip(d_list,
            merged_data_mapfiles)]
        image1_basenames_mapfiles = self.s.run(actions)

        self.log.info('FFTing model image (facet model #1)...')
        self.log.debug('Merging unaverged chunks...')
        merged_unavg_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            unavg_files, _ = read_mapfile(chunk_data_mapfiles[i])
            merged_unavg_data_mapfiles.append(self.write_mapfile([merge_chunks(unavg_files,
            prefix=None, clobber=True)], direction=d_list[i],
            host_list=d_hosts[i]))
        actions = [FFT(self.parset, dm, mm, p['imager1'], prefix='fft1',
            direction=d, index=1) for d, dm, mm in zip(d_list,
            merged_unavg_data_mapfiles, image1_basenames_mapfiles)]
        self.s.run(actions)

        chunk_parmdb_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_parmdb_mapfiles.append(self.write_mapfile(
                [chunk.parmdb_phaseonly2 for chunk in chunks],
                prefix='chunk_parmdb_phaseonly2', direction=d_list[i],
                host_list=d_hosts[i]))

        self.log.info('Solving for phase solutions and applying them (#2)...')
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseonly2'].copy()
            p_d['timestep'] = d.solint_p
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=None,
            parmdb_datamap=pm, prefix='facet_phaseonly', direction=d, index=1)
            for d, dm, pd, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_parmdb_mapfiles)]
        self.s.run(actions)

        self.log.info('Imaging (facet image #2)...')
        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg2'], prefix='facet',
            direction=d, index=2) for d, m in zip(d_list, chunk_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Merging averged chunks...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            avg_files, _ = read_mapfile(avg_data_mapfiles[i])
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(avg_files,
            prefix=None, clobber=True)], direction=d_list[i],
            host_list=d_hosts[i]))

        self.log.debug('Imaging...')
        actions = [MakeImageIterate(self.parset, m, p['imager2'], prefix='facet_selfcal2',
            direction=d) for d, m in zip(d_list, merged_data_mapfiles)]
        image2_basenames_mapfiles = self.s.run(actions)

        self.log.info('FFTing model image (facet model #2)...')
        actions = [FFT(self.parset, dm, mm, p['imager2'], prefix='fft2',
            direction=d, index=2) for d, dm, mm in zip(d_list,
            merged_unavg_data_mapfiles, image2_basenames_mapfiles)]
        self.s.run(actions)

        chunk_parmdb_phaseamp_phase1_mapfiles = []
        chunk_parmdb_phaseamp_amp1_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_parmdb_phaseamp_phase1_mapfiles.append(self.write_mapfile(
                [chunk.parmdb_phaseamp_phase1 for chunk in chunks],
                prefix='chunk_parmdb_phase1', direction=d_list[i],
                host_list=d_hosts[i]))
            chunk_parmdb_phaseamp_amp1_mapfiles.append(self.write_mapfile(
                [chunk.parmdb_phaseamp_amp1 for chunk in chunks],
                prefix='chunk_parmdb_amp1', direction=d_list[i],
                host_list=d_hosts[i]))

        self.log.info('Solving for amplitude solutions and applying them (#1)...')
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseamp1_phaseonly'].copy()
            p_d['timestep'] = d.solint_p
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=None,
            parmdb_datamap=pm, prefix='facet_phaseonly', direction=d, index=2)
            for d, dm, pd, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_parmdb_phaseamp_phase1_mapfiles)]
        self.s.run(actions)
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseamp1_amponly'].copy()
            p_d['timestep'] = d.solint_a
            p_d['chunksize'] = d.solint_a
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=None,
            parmdb_datamap=pm, prefix='facet_amponly', direction=d, index=2)
            for d, dm, pd, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_parmdb_phaseamp_amp1_mapfiles)]
        self.s.run(actions)

        self.log.debug('Merging chunk instrument parmdbs...')
        merged_parmdb_phaseamp_amp1_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            merged_parmdb_phaseamp_amp1_mapfiles.append(self.write_mapfile(
                [merge_chunk_parmdbs([chunk.parmdb_phaseamp_amp1 for chunk
                in chunks], prefix='merged_chunks_amps1')], prefix='merged_amps1',
                direction=d_list[i], host_list=d_hosts[i]))

        self.log.info('Smoothing amplitude solutions...')
        actions = [Smooth(self.parset, dm, p['smooth_amp1'], pm,
            prefix='facet_amp', direction=d, index=2)
            for d, dm, pm in zip(d_list, facet_data_mapfiles,
            merged_parmdb_phaseamp_amp1_mapfiles)]
        self.s.run(actions)

        self.log.info('Applying amplitude solutions...')
        actions = [Apply(self.parset, dm, p['apply_amp1'],
            pm, prefix='facet_amp', direction=d, index=2) for d, dm, pm in
            zip(d_list, chunk_data_mapfiles,
            merged_parmdb_phaseamp_amp1_mapfiles)]
        self.s.run(actions)

        self.log.info('Imaging (facet image #3)...')
        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg3'], prefix='facet',
            direction=d, index=3) for d, m in zip(d_list, chunk_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            avg_files, _ = read_mapfile(avg_data_mapfiles[i])
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(avg_files,
            prefix=None, clobber=True)], direction=d_list[i],
            host_list=d_hosts[i]))

        self.log.debug('Imaging...')
        actions = [MakeImageIterate(self.parset, m, p['imager3'], prefix='facet_selfcal3',
            direction=d) for d, m in zip(d_list, merged_data_mapfiles)]
        image3_basenames_mapfiles = self.s.run(actions)

        self.log.info('FFTing model image (facet model #3)...')
        actions = [FFT(self.parset, dm, mm, p['imager3'], prefix='fft3',
            direction=d, index=3) for d, dm, mm in zip(d_list,
            merged_unavg_data_mapfiles, image3_basenames_mapfiles)]
        self.s.run(actions)

        self.log.info('Solving for amplitude solutions (#2)...')
        chunk_parmdb_phaseamp_phase2_mapfiles = []
        chunk_parmdb_phaseamp_amp2_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_parmdb_phaseamp_phase2_mapfiles.append(self.write_mapfile(
                [chunk.parmdb_phaseamp_phase2 for chunk in chunks],
                prefix='chunk_parmdb_phase2', direction=d_list[i],
                host_list=d_hosts[i]))
            chunk_parmdb_phaseamp_amp2_mapfiles.append(self.write_mapfile(
                [chunk.parmdb_phaseamp_amp2 for chunk in chunks],
                prefix='chunk_parmdb_amp2', direction=d_list[i],
                host_list=d_hosts[i]))

        p_list = []
        for d in d_list:
            p_d = p['solve_phaseamp2_phaseonly'].copy()
            p_d['timestep'] = d.solint_p
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=None,
            parmdb_datamap=pm, prefix='facet_phaseonly', direction=d, index=3)
            for d, dm, pd, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_parmdb_phaseamp_phase2_mapfiles)]
        self.s.run(actions)
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseamp2_amponly'].copy()
            p_d['timestep'] = d.solint_a
            p_d['chunksize'] = d.solint_a
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=None,
            parmdb_datamap=pm, prefix='facet_amponly', direction=d, index=3)
            for d, dm, pd, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_parmdb_phaseamp_amp2_mapfiles)]
        self.s.run(actions)

        self.log.debug('Merging chunk instrument parmdbs...')
        merged_parmdb_phaseamp_amp2_mapfiles = []
        merged_parmdb_phaseamp_phase2_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            merged_parmdb_phaseamp_amp2_mapfiles.append(self.write_mapfile(
            	[merge_chunk_parmdbs([chunk.parmdb_phaseamp_amp2 for
            	chunk in chunks], prefix='merged_amps2')],
            	prefix='merged_amps2', direction=d_list[i], host_list=d_hosts[i]))
            merged_parmdb_phaseamp_phase2_mapfiles.append(self.write_mapfile(
            	[merge_chunk_parmdbs([chunk.parmdb_phaseamp_phase2 for
            	chunk in chunks], prefix='merged_phases2')],
            	prefix='merged_phases2', direction=d_list[i], host_list=d_hosts[i]))

        self.log.info('Smoothing amplitude solutions...')
        actions = [Smooth(self.parset, dm, p['smooth_amp2'], pm,
            prefix='facet_amp', direction=d, index=3)
            for d, dm, pm in zip(d_list, facet_data_mapfiles,
            merged_parmdb_phaseamp_amp2_mapfiles)]
        self.s.run(actions)

        self.log.info('Applying amplitude solutions...')
        actions = [Apply(self.parset, dm, p['apply_amp3'],
            pm, prefix='facet_amp', direction=d, index=4) for d, dm, pm in
            zip(d_list, chunk_data_mapfiles,
            merged_parmdb_phaseamp_amp2_mapfiles)]
        self.s.run(actions)

        self.log.info('Imaging (facet image #4)...')
        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg4'], prefix='facet',
            direction=d, index=4) for d, m in zip(d_list, chunk_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            avg_files, _ = read_mapfile(avg_data_mapfiles[i])
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(avg_files,
            prefix=None, clobber=True)], direction=d_list[i],
            host_list=d_hosts[i]))

        self.log.debug('Imaging...')
        actions = [MakeImageIterate(self.parset, m, p['imager4'], prefix='facet_selfcal4',
            direction=d) for d, m in zip(d_list, merged_data_mapfiles)]
        image4_basenames_mapfiles = self.s.run(actions)

        # Check if image rms / ration of max/min is still improving. If so,
        # continue last selfcal step. If not, stop sefcal

        # Check images if interactive=True
        # If not OK: stop and reset state for this direction
        # If OK, continue
        if self.parset['interactive']:
            import pyrap.images as pim
            for i, d in enumerate(d_list):
                self.log.info('Showing selfcal images (initial, 2nd phase-only, '
                    '2nd phase+amp) for direction {0}...'.format(d.name))
                image1, _ = read_mapfile(image0_basenames_mapfiles[i])
                image2, _ = read_mapfile(image2_basenames_mapfiles[i])
                image3, _ = read_mapfile(image4_basenames_mapfiles[i])
                im = pim.image([image1[0], image2[0], image3[0]])
                im.view()
                prompt = "Continue processing (y/n)? "
                answ = raw_input(prompt)
                while answ.lower() not in  ['y', 'n', 'yes', 'no']:
                    answ = raw_input(prompt)
                if answ.lower() in ['n', 'no']:
                    self.log.info('Exiting...')
                    sys.exit()

        self.log.info('Merging final instrument parmdbs...')
        merged_parmdb_final_mapfiles = []
        for i, d in enumerate(d_list):
            phases2_final, _ = read_mapfile(merged_parmdb_phaseamp_phase2_mapfiles[i])
            smoothed_amps2_final, _ = read_mapfile(merged_parmdb_phaseamp_amp2_mapfiles[i])
            merged_file = merge_parmdbs(phases2_final[0], smoothed_amps2_final[0])
            merged_parmdb_final_mapfiles.append(self.write_mapfile([merged_file],
                prefix='merged_amps_phases_final', direction=d, host_list=d_hosts[i]))

        # Save files to the direction objects
        for d, m in zip(d_list, merged_parmdb_final_mapfiles):
            f, _ = read_mapfile(m)
            d.dirdepparmdb = f[0]

        # Save state
        self.set_completed(d_list)


class FacetImage(Operation):
    """
    Operation to image the full facet
    """
    def __init__(self, parset, bands, direction=None, reset=False,
        name='FacetImage'):
        super(FacetImage, self).__init__(parset, bands, direction=direction,
            reset=reset, name=name)

        # Set up scheduler (runs at most num_nodes directions in parallel)
        num_nodes = len(self.parset['cluster_specific']['node_list'])
        self.s = Scheduler(max_threads=num_nodes, name=self.name,
            op_parset=self.parset)


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import Average, Concatenate, ChgCentre
        from factor.actions.calibrations import Apply
        from factor.actions.images import MakeImageIterate
        from factor.actions.models import FFT
        from factor.lib.operation_lib import copy_column, merge_chunks
        from factor.operations.hardcoded_param import facet_image as p
        from factor.lib.datamap_lib import read_mapfile

        bands = self.bands
        d_list = self.direction
        if type(d_list) is not list:
            d_list = [d_list]

        # Check state
        if self.check_completed(d_list):
            return

        # Set image sizes
        for d in d_list:
            cell = float(p['imager']['cell'].split('arcsec')[0]) # arcsec per pixel
            imsize = d.width * 1.1 * 3600.0 / cell # pixels
            if imsize < 512:
                imsize = 512
            d.imsize = imsize

        # Divide up the nodes among the directions
        node_list = self.parset['cluster_specific']['node_list']
        if len(d_list) >= len(node_list):
            for i in range(len(d_list)-len(node_list)):
                node_list.append(node_list[i])
            d_hosts = [[n] for n in node_list]
        else:
            parts = len(d_list)
            d_hosts = [node_list[i*len(node_list)//parts:
                (i+1)*len(node_list)//parts] for i in range(parts)]

        # Make initial data maps for the phase-shifted datasets and their dir-dep
        # instrument parmdbs
        shifted_all_data_mapfiles = []
        dir_dep_parmdbs_mapfiles = []
        for d, h in zip(d_list, d_hosts):
            shifted_all_data_mapfiles.append(self.write_mapfile(d.shifted_all_data_files,
                prefix='shifted', direction=d, host_list=h))
            dir_dep_parmdbs_mapfiles.append(self.write_mapfile([d.dirdepparmdb]*
                len(bands), prefix='dir_dep_parmdbs', direction=d, host_list=h))

        self.log.info('Applying direction-dependent calibration...')
        actions = [Apply(self.parset, dm, p['apply_dirdep'],
            pm, prefix='facet_dirdep', direction=d) for d, dm, pm in zip(d_list,
            shifted_all_data_mapfiles, dir_dep_parmdbs_mapfiles)]
        self.s.run(actions)

        self.log.info('Averaging...')
        actions = [Average(self.parset, m, p['avg'], prefix='facet',
            direction=d) for d, m in zip(d_list, shifted_all_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.info('Imaging...')
        if self.parset['use_chgcentre']:
            self.log.debug('Changing center to zenith...')
            actions = [ChgCentre(self.parset, dm, {},
                prefix='facet', direction=d) for d, dm in
                zip(d_list, avg_data_mapfiles)]
            chgcentre_data_mapfiles = self.s.run(actions)
            input_to_imager_mapfiles = chgcentre_data_mapfiles
        else:
            input_to_imager_mapfiles = avg_data_mapfiles

        self.log.debug('Merging files...')
        merged_data_mapfiles = []
        for i, mf in enumerate(input_to_imager_mapfiles):
            in_files, _ = read_mapfile(mf)
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(in_files,
            prefix=None, clobber=True)], prefix='imager', direction=d_list[i],
            host_list=d_hosts[i]))

        actions = [MakeImageIterate(self.parset, m, p['imager'], prefix='facet_image',
            direction=d) for d, m in zip(d_list, merged_data_mapfiles)]
        image_basenames_mapfiles = self.s.run(actions)

        self.log.info('FFTing model image...')
        merged_data_mapfiles = []
        for i, mf in enumerate(shifted_all_data_mapfiles):
            in_files, _ = read_mapfile(mf)
            merged_data_mapfiles.append(self.write_mapfile([merge_chunks(in_files,
            prefix=None, clobber=True)], direction=d_list[i],
            host_list=d_hosts[i]))
        actions = [FFT(self.parset, dm, mm, p['imager'], prefix='fft',
            direction=d) for d, dm, mm in zip(d_list,
            merged_data_mapfiles, image_basenames_mapfiles)]
        self.s.run(actions)

        # Save image basenames to the direction objects
        for d, mf in zip(d_list, image_basenames_mapfiles):
            file, _ = read_mapfile(mf)
            d.skymodel_dirdep = file[0]

        # Save state
        self.set_completed(d_list)


class FacetSub(Operation):
    """
    Operation to subtract final facet sky model from facet visibilties
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetSub, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetSub')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.calibrations import Subtract
        from factor.actions.visibilities import PhaseShift, Average, ChgCentre
        from factor.actions.images import MakeImage
        from factor.operations.hardcoded_param import facet_sub as p
        from factor.lib.datamap_lib import read_mapfile
        from factor.lib.operation_lib import copy_column, verify_subtract

        bands = self.bands
        d_list = self.direction
        if type(d_list) is not list:
            d_list = [d_list]

        # Check state
        if self.check_completed(d_list):
            return

        # Divide up the nodes among the directions
        node_list = self.parset['cluster_specific']['node_list']
        if len(d_list) >= len(node_list):
            for i in range(len(d_list)-len(node_list)):
                node_list.append(node_list[i])
            d_hosts = [[n] for n in node_list]
        else:
            parts = len(d_list)
            d_hosts = [node_list[i*len(node_list)//parts:
                (i+1)*len(node_list)//parts] for i in range(parts)]

        # Make initial data maps for the empty datasets, their dir-dep
        # instrument parmdbs, and their dir-dep sky models
        shifted_all_data_mapfiles = []
        shifted_sub_data_mapfiles = []
        dir_dep_parmdbs_mapfiles = []
        for d, h in zip(d_list, d_hosts):
            shifted_all_data_mapfiles.append(self.write_mapfile(d.shifted_all_data_files,
            	prefix='shifted_all', direction=d, host_list=h))
            shifted_sub_data_mapfiles.append(self.write_mapfile(d.shifted_sub_data_files,
            	prefix='shifted_sub', direction=d, host_list=h))
            dir_dep_parmdbs_mapfiles.append(self.write_mapfile([d.dirdepparmdb]*
                len(bands), prefix='dir_dep_parmdbs', direction=d, host_list=h))

        self.log.info('Subtracting sources...')
        actions = [Subtract(self.parset, dm, p['subtract'],
            model_datamap=None, parmdb_datamap=pd, prefix='facet_dirdep',
            direction=d) for d, dm, pd in
            zip(d_list, shifted_all_data_mapfiles, dir_dep_parmdbs_mapfiles)]
        self.s.run(actions)

        # apply direction-independent calibration
        self.log.info('Applying direction-independent calibration...')
        actions = [Apply(self.parset, dm, p['apply_pre'],
            pm, prefix='facet_dirindep', direction=d, index=1) for d, dm, pm in zip(d_list,
            shifted_sub_data_mapfiles, dir_indep_parmdbs_mapfiles)]
        self.s.run(actions)
        actions = [Apply(self.parset, dm, p['apply_post'],
            pm, prefix='facet_dirindep', direction=d, index=2) for d, dm, pm in zip(d_list,
            shifted_all_data_mapfiles, dir_indep_parmdbs_mapfiles)]
        self.s.run(actions)

        self.log.info('Phase shifting back to field center...')
        ra = bands[0].ra
        dec = bands[0].dec
        actions = [PhaseShift(self.parset, dm, p['shift'], prefix='facet',
            direction=d, ra=ra, dec=dec, index=1) for d, dm in zip(d_list,
            shifted_sub_data_mapfiles)]
        unshifted_pre_data_mapfiles = self.s.run(actions)
        actions = [PhaseShift(self.parset, dm, p['shift'], prefix='facet',
            direction=d, ra=ra, dec=dec, index=2) for d, dm in zip(d_list,
            shifted_all_data_mapfiles)]
        unshifted_post_data_mapfiles = self.s.run(actions)

        self.log.info('Averaging...')
        actions = [Average(self.parset, dm, p['avg'], prefix='facet', index=1)
            for d, dm in zip(d_list, unshifted_pre_data_mapfiles)]
        avg_pre_data_mapfiles = self.s.run(actions)
        actions = [Average(self.parset, dm, p['avg'], prefix='facet', index=2)
            for d, dm in zip(d_list, unshifted_post_data_mapfiles)]
        avg_post_data_mapfiles = self.s.run(actions)

        self.log.info('Imaging...')
        if self.parset['use_chgcentre']:
            self.log.debug('Changing center to zenith...')
            actions = [ChgCentre(self.parset, dm, {},
                prefix='resid', direction=d, index=1) for d, dm in
                zip(d_list, avg_data_mapfiles)]
            chgcentre_pre_data_mapfiles = self.s.run(actions)
            input_to_imager_pre_mapfiles = chgcentre_pre_data_mapfiles
            actions = [ChgCentre(self.parset, dm, {},
                prefix='resid', direction=d, index=2) for d, dm in
                zip(d_list, avg_data_mapfiles)]
            chgcentre_post_data_mapfiles = self.s.run(actions)
            input_to_imager_post_mapfiles = chgcentre_post_data_mapfiles
        else:
            input_to_imager_pre_mapfiles = avg_pre_data_mapfiles
            input_to_imager_post_mapfiles = avg_post_data_mapfiles
        actions = [MakeImage(self.parset, dm, p['imager'], prefix='field_image',
            direction=d, index=1) for d, dm in zip(d_list, input_to_imager_pre_mapfiles)]
        image_pre_basenames_mapfiles = self.s.run(actions)
        actions = [MakeImage(self.parset, dm, p['imager'], prefix='field_image',
            direction=d, index=2) for d, dm in zip(d_list, input_to_imager_post_mapfiles)]
        image_post_basenames_mapfiles = self.s.run(actions)

        self.log.info('Checking residual images...')
        for d, impre_mapfile, impost_mapfile in zip(d_list, image_pre_basenames_mapfiles,
            image_post_basenames_mapfiles):
            image_pre_files, _ = read_mapfile(impre_mapfile)
            image_post_files, _ = read_mapfile(impost_mapfile)
            res_val = 0.5
            # Check only lowest-frequency images for now
            d.selfcal_ok = verify_subtract(image_pre_files[0], image_post_files[0],
                res_val, self.parset['imager'])

        # Save state
        self.set_completed(d_list)


class FacetAddAllFinal(Operation):
    """
    Operation to add all sources in the facet to data (final)
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetAddAllFinal, self).__init__(parset, bands, direction=direction,
            reset=reset, name=name)


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import PhaseShift
        from factor.actions.calibrations import Add
        from factor.actions.models import FFT
        from factor.operations.hardcoded_param import facet_add_all_final as p
        from factor.lib.datamap_lib import read_mapfile

        bands = self.bands
        d = self.direction

        # Check state
        if self.check_completed(d):
            return

        # Make initial data maps for the empty datasets, their dir-indep
        # instrument parmdbs, and their dir-indep sky models
        subtracted_all_mapfile = self.write_mapfile([band.file for band in bands],
        	prefix='subtracted_all', direction=d)
        dir_dep_parmdbs_mapfile = self.write_mapfile([d.dirdepparmdb]*len(bands),
            prefix='dir_dep_parmdbs', direction=d)
        dir_dep_models_mapfile = self.write_mapfile([d.skymodel_dirdep]*len(bands),
            prefix='dir_dep_skymodels', direction=d)

        self.log.debug('FFTing model image (facet model final)...')
        action = FFT(self.parset, subtracted_all_mapfile,
            dir_dep_models_mapfile, p['fft'], prefix='fft', direction=d)
        self.s.run(action)

        self.log.info('Adding sources for this direction...')
        action = Add(self.parset, subtracted_all_mapfile, p['add'],
            dir_indep_cal_skymodels_mapfile, dir_indep_parmdbs_mapfile,
            prefix='facet_dirdep', direction=d)
        self.s.run(action)

        self.log.info('Phase shifting DATA...')
        action = PhaseShift(self.parset, subtracted_all_mapfile, p['shift'],
            prefix='facet', direction=d)
        shifted_data_mapfile = self.s.run(action)

        # Save files to the direction object
        d.shifted_data_files, _ = read_mapfile(shifted_data_mapfile)

        # Save state
        self.set_completed(d)


def FacetImageFinal(FacetImage):
    """
    Operation to make final facet image
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetImageFinal, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetImageFinal')

