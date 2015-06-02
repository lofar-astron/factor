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
FacetCheck : Operation
    Subtracts all facet sources from data and images residuals
FacetSub : Operation
    Subtracts all facet sources from data
FacetAddAllFinal : Operation
    Adds all facet sources from final model
FacetImageFinal : Operation
    Images the entire facet

"""
import os
from factor.lib.operation import Operation
from factor.lib.scheduler_mp import Scheduler


class FacetAdd(Operation):
    """
    Operation to add calibrator source to data
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetAdd, self).__init__(parset, bands, direction,
            name='FacetAdd')

        # Define parameters needed for this operation
        skymodels = [band.skymodel_dirindep for band in self.bands]
        self.parms_dict = {'input_dir': parset['dir_ms'],
                           'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'dir_indep_parmdb_name': parset['parmdb_name'],
                           'skymodels': skymodels,
                           'facet_ra': self.direction.ra,
                           'facet_dec': self.direction.dec,
                           'cal_radius_deg': self.direction.cal_radius_deg,
                           'facet_state_file': self.direction.save_file,
                           'hosts': self.node_list}


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to direction object
        self.direction.shifted_all_bands_datamap = os.path.join(self.mapfile_dir,
            'shifted_all_bands.datamap')
        self.direction.shifted_cal_bands_datamap = os.path.join(self.mapfile_dir,
            'shifted_all_bands.datamap')
        self.direction.shifted_empty_bands_datamap = os.path.join(self.mapfile_dir,
            'shifted_all_bands.datamap')
        self.direction.dir_indep_parmdbs_datamap = os.path.join(self.mapfile_dir,
            'dir_indep_parmdbs.datamap')


class FacetSetup(Operation):
    """
    Operation to set up data for selfcal
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetSetup, self).__init__(parset, bands, direction,
            name='FacetSetup')

        # Define parameters needed for this operation
        self.parms_dict = {'input_dir': parset['dir_ms'],
                           'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'dir_indep_parmdbs_datamap': self.direction.dir_indep_parmdbs_datamap,
                           'hosts': self.direction.hosts}


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamap to direction object
        self.direction.shifted_cal_concat_datamap = os.path.join(self.mapfile_dir,
            'shifted_cal_concat_bands.datamap')


class FacetSelfcal(Operation):
    """
    Operation to selfcal one or more directions
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetSelfcal, self).__init__(parset, bands, direction,
            name='FacetSelfcal')

        # Define parameters needed for this operation
        self.parms_dict = {'input_dir': parset['dir_ms'],
                           'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'shifted_cal_concat_datamap': self.direction.shifted_cal_concat_datamap,
                           'hosts': self.direction.hosts}


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamap to direction object
        self.direction.shifted_cal_concat_datamap = os.path.join(self.mapfile_dir,
            'shifted_cal_concat_bands.datamap')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import Average, Concatenate, PhaseShift
        from factor.actions.calibrations import Apply, Solve
        from factor.actions.images import MakeImage, MakeMask, image_with_mask
        from factor.actions.models import FFT
        from factor.actions.solutions import Smooth, ResetPhases
        from factor.lib.operation_lib import copy_column, make_chunks, \
            merge_chunks, merge_parmdbs, merge_chunk_parmdbs, check_selfcal
        from factor.operations.hardcoded_param import facet_selfcal as p
        from factor.lib.datamap_lib import read_mapfile
        import numpy as np

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
        image0_basenames_mapfiles = image_with_mask(self, p['imager0'],
            'facet_selfcal0', avg_data_mapfiles, directions=d_list)

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
        image1_basenames_mapfiles = image_with_mask(self, p['imager1'],
            'facet_selfcal1', merged_data_mapfiles, directions=d_list)

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
        image2_basenames_mapfiles = image_with_mask(self, p['imager2'],
            'facet_selfcal2', merged_data_mapfiles, directions=d_list)

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
        image_final_basenames_mapfiles = image_with_mask(self, p['imager3'],
            'facet_selfcal3', merged_data_mapfiles, directions=d_list)

        # Loop over final calibration as long as there is improvement
        index = 3
        for d in d_list:
            d.improving = True
        while np.any([d.improving for d in d_list]):
            self.log.info('FFTing model image (facet model #3)...')
            actions = [FFT(self.parset, dm, mm, p['imager3'], prefix='fft3',
                direction=d, index=index) for d, dm, mm in zip(d_list,
                merged_unavg_data_mapfiles, image_final_basenames_mapfiles)
                if d.improving]
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
                parmdb_datamap=pm, prefix='facet_phaseonly', direction=d, index=index)
                for d, dm, pd, pm in zip(d_list, chunk_data_mapfiles, p_list,
                chunk_parmdb_phaseamp_phase2_mapfiles) if d.improving]
            self.s.run(actions)
            p_list = []
            for d in d_list:
                p_d = p['solve_phaseamp2_amponly'].copy()
                p_d['timestep'] = d.solint_a
                p_d['chunksize'] = d.solint_a
                p_list.append(p_d)
            actions = [Solve(self.parset, dm, pd, model_datamap=None,
                parmdb_datamap=pm, prefix='facet_amponly', direction=d, index=index)
                for d, dm, pd, pm in zip(d_list, chunk_data_mapfiles, p_list,
                chunk_parmdb_phaseamp_amp2_mapfiles) if d.improving]
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
                prefix='facet_amp', direction=d, index=index)
                for d, dm, pm in zip(d_list, facet_data_mapfiles,
                merged_parmdb_phaseamp_amp2_mapfiles) if d.improving]
            self.s.run(actions)

            index += 1

            self.log.info('Applying amplitude solutions...')
            actions = [Apply(self.parset, dm, p['apply_amp3'],
                pm, prefix='facet_amp', direction=d, index=index) for d, dm, pm in
                zip(d_list, chunk_data_mapfiles,
                merged_parmdb_phaseamp_amp2_mapfiles) if d.improving]
            self.s.run(actions)

            self.log.info('Imaging (facet image #4)...')
            self.log.debug('Averaging in preparation for imaging...')
            actions = [Average(self.parset, m, p['avg4'], prefix='facet',
                direction=d, index=index) for d, m in zip(d_list, chunk_data_mapfiles)
                if d.improving]
            avg_data_mapfiles = self.s.run(actions)

            self.log.debug('Merging chunks...')
            merged_data_mapfiles = []
            for i, chunks in enumerate(chunks_list):
                avg_files, _ = read_mapfile(avg_data_mapfiles[i])
                merged_data_mapfiles.append(self.write_mapfile([merge_chunks(avg_files,
                    prefix=None, clobber=True)], direction=d_list[i],
                    host_list=d_hosts[i]))

            self.log.debug('Imaging...')
            image_final_basenames_mapfiles_prev = image_final_basenames_mapfiles[:]
            d_list_impr = []
            merged_data_mapfiles_impr = []
            for d, m in zip(d_list, merged_data_mapfiles):
                if d.improving:
                    d_list_impr.append(d)
                    merged_data_mapfiles_impr.append(m)
            if len(d_list_impr) > 0:
                image_final_basenames_mapfiles = image_with_mask(self, p['imager4'],
                    'facet_selfcal{0}'.format(index), merged_data_mapfiles_impr,
                    directions=d_list_impr)

            # Check if image rms / ratio of max to min are still improving. If so,
            # continue last selfcal step. If not, stop sefcal
            for d, pm, fm in zip(d_list, image_final_basenames_mapfiles_prev,
                image_final_basenames_mapfiles):
                d.loop_amp_selfcal = False # not yet implemented!
                if d.improving:
                    if d.loop_amp_selfcal:
                        prev_image, _ = read_mapfile(pm)
                        new_image, _ = read_mapfile(fm)
                        d.improving = check_selfcal(prev_image[0], new_image[0],
                            self.p['check_selfcal']['rms_threshold'],
                            self.p['check_selfcal']['ratio_threshold'])
                    else:
                        d.improving = False

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
        self.s = Scheduler(max_procs=num_nodes, name=self.name,
            op_parset=self.parset)


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import Average, Concatenate, ChgCentre
        from factor.actions.calibrations import Apply
        from factor.actions.images import image_with_mask
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

        image_basenames_mapfiles = image_with_mask(self, p['imager'],
            'facet_image', merged_data_mapfiles, directions=d_list)

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


class FacetCheck(Operation):
    """
    Operation to subtract final facet sky model from facet visibilties
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetCheck, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetCheck')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.calibrations import Subtract, Apply
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
        dir_indep_parmdbs_mapfiles = []
        dir_dep_parmdbs_mapfiles = []
        for d, h in zip(d_list, d_hosts):
            shifted_all_data_mapfiles.append(self.write_mapfile(d.shifted_all_data_files,
            	prefix='shifted_all', direction=d, host_list=h))
            shifted_sub_data_mapfiles.append(self.write_mapfile(d.shifted_sub_data_files,
            	prefix='shifted_sub', direction=d, host_list=h))
            dir_dep_parmdbs_mapfiles.append(self.write_mapfile([d.dirdepparmdb]*
                len(bands), prefix='dir_dep_parmdbs', direction=d, host_list=h))
            dir_indep_parmdbs_mapfiles.append(self.write_mapfile([band.dirindparmdb
                for band in bands], prefix='dir_indep_parmdbs', direction=d, host_list=h))

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
            d.max_residual_val = 0.25
            image_pre_files, _ = read_mapfile(impre_mapfile)
            image_post_files, _ = read_mapfile(impost_mapfile)
            # Check only lowest-frequency images for now
            d.selfcal_ok, maxval, maxvalpre = verify_subtract(image_pre_files[0], image_post_files[0],
                d.max_residual_val, self.parset['imager'])
            self.log.info('Current residual: {0} Jy/beam; previous residual: '
                '{1} Jy/beam'.format(maxval, maxvalpre))

        # Save state
        self.set_completed(d_list)


class FacetSub(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetSub, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetSub')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import PhaseShift
        from factor.actions.calibrations import Subtract
        from factor.lib.datamap_lib import read_mapfile
        from factor.lib.operation_lib import copy_column
        from factor.operations.hardcoded_param import field_sub as p

        bands = self.bands
        d = self.direction

        # Check state
        if self.check_completed(d):
            return

        # Make initial data maps for the input datasets and their
        # instrument parmdbs.
        shifted_all_data_mapfile = self.write_mapfile(d.shifted_all_data_files,
            prefix='shifted', direction=d)
        dir_dep_parmdbs_mapfile = self.write_mapfile([d.dirdepparmdb]*len(bands),
            prefix='dir_dep_parmdbs', direction=d)
        orig_data_mapfile = self.write_mapfile([band.file for band in bands],
        	prefix='input_data')

        self.log.info('Phase shifting model back to field center...')
        ra = bands[0].ra
        dec = bands[0].dec
        action = PhaseShift(self.parset, shifted_all_data_mapfile, p['shift'],
            prefix='facet', direction=d, ra=ra, dec=dec)
        unshifted_model_data_mapfile = self.s.run(action)

        self.log.info('Copying model to bands...')
        ms_to_files = [band.file for band in bands]
        ms_from_files, _ = read_mapfile(unshifted_model_data_mapfile)
        for ms_to, ms_from in zip(ms_to_files, ms_from_files):
            copy_column(ms_to, p['copy']['incol'], p['copy']['outcol'], ms_from)

        self.log.info('Subtracting final model for this direction...')
        action = Subtract(self.parset, orig_data_mapfile, p['subtract'], None,
        	dir_dep_parmdbs_mapfile, prefix='field_dirdep', direction=d)
        self.s.run(action)

        # Save state
        self.set_completed(d)


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

