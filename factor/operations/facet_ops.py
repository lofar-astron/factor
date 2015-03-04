"""
Module that holds all facet operations

Classes
-------
FacetAddCal : Operation
    Adds calibrator source to data
FacetSetup : Operation
    Sets up the data for selfcal
FacetSelfcal : Operation
    Runs selfcal cycle
FacetAddAll : Operation
    Adds all facet sources to data
FacetSetupFull : Operation
    Sets up the data for full imaging
FacetImage : Operation
    Images the entire facet
FacetSubAll : Operation
    Subtracts all facet sources from data
FacetAddAllFinal : Operation
    Adds all facet sources to data
FacetImageFinal : Operation
    Images the entire facet

"""
import os
from factor.lib.operation import Operation


class FacetAddCal(Operation):
    """
    Operation to add calibrator source to data
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetAddCal, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetAddCal')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import PhaseShift
        from factor.actions.models import MakeFacetSkymodel
        from factor.actions.calibrations import Add
        from factor.operations.hardcoded_param import facet_add_cal as p
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        d = self.direction
        bands = self.bands

        if os.path.exists(self.statebasename+'.done'):
            shifted_data_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/FacetAddCal/PhaseShift/facet_output_{0}.datamap'.
                format(d.name))
            files = read_mapfile(shifted_data_mapfile)
            for band, f in zip(bands, files):
                band.shifted_data_file = f
            return

        # Make initial data maps for the empty datasets, their dir-indep
        # instrument parmdbs, and their dir-indep sky models
        subtracted_all_mapfile = write_mapfile([band.file for band in bands],
        	self.name, prefix='subtracted_all',
        	working_dir=self.parset['dir_working'])
        dir_indep_parmdbs_mapfile = write_mapfile([band.dirindparmdb for band
        	in bands], self.name, prefix='dir_indep_parmdbs',
        	working_dir=self.parset['dir_working'])
        dir_indep_skymodels_mapfile = write_mapfile([band.skymodel_dirindep
        	for band in bands], self.name, prefix='dir_indep_skymodels',
        	working_dir=self.parset['dir_working'])

        # Add calibrators from the dir-indep sky model for this direction to the
        # visibilities
        self.log.info('Selecting sources for this direction...')
        action = MakeFacetSkymodel(self.parset, dir_indep_skymodels_mapfile,
            p['select'], d, prefix='cal', cal_only=True)
        dir_indep_cal_skymodels_mapfile = action.run()

        self.log.info('Adding sources for this direction...')
        self.parset['use_ftw'] = False
        action = Add(self.parset, subtracted_all_mapfile, p['add'],
            dir_indep_cal_skymodels_mapfile, dir_indep_parmdbs_mapfile,
            prefix='facet_dirindep', direction=d)
        action.run()

        # Phase shift to facet center
        self.log.info('Phase shifting DATA...')
        action = PhaseShift(self.parset, subtracted_all_mapfile, p['shift'],
            prefix='facet', direction=d)
        shifted_data_mapfile = action.run()

        # Save files to the band objects
        files = read_mapfile(shifted_data_mapfile)
        for band, f in zip(bands, files):
            band.shifted_data_file = f


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
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        d_list = self.direction
        if type(d_list) is not list:
            d_list = [d_list]
        bands = self.bands

        # Check state for each direction
        all_done = False
        for i, d in enumerate(d_list):
            if os.path.exists(self.statebasename[i]+'.done'):
                concat_data_mapfile = os.path.join(self.parset['dir_working'],
                    'datamaps/FacetSetup/Concatenate/facet_bands_output_{0}.datamap'.
                    format(d.name))
                file = read_mapfile(concat_data_mapfile)[0]
                d.concat_file = file
                all_done = True
        if all_done:
            return

        # Make initial data maps for the phase-shifted datasets and their dir-indep
        # instrument parmdbs
        shifted_data_mapfiles = []
        dir_indep_parmdbs_mapfiles = []
        for d in d_list:
            shifted_data_mapfiles.append(write_mapfile([band.shifted_data_file
                for band in bands], self.name, prefix='shifted',
                working_dir=self.parset['dir_working']))
            dir_indep_parmdbs_mapfiles.append(write_mapfile([band.dirindparmdb
                for band in bands], self.name, prefix='dir_indep_parmdbs',
                working_dir=self.parset['dir_working']))

        # apply direction-independent calibration
        self.log.info('Applying direction-independent calibration...')
        actions = [Apply(self.parset, dm, p['apply'],
            pm, prefix='facet_dirindep', direction=d) for d, dm, pm in zip(d_list,
            shifted_data_mapfiles, dir_indep_parmdbs_mapfiles)]
        self.s.run(actions)

        # average to 1 channel per band. Do this twice, once for DATA and once
        # for CORRECTED_DATA
        self.log.info('Averaging DATA...')
        actions = [Average(self.parset, m, p['avg1'], prefix='facet',
            direction=d, index=1) for d, m in zip(d_list, shifted_data_mapfiles)]
        avg1_data_mapfiles = self.s.run(actions)
        self.log.info('Averaging CORRECTED_DATA...')
        actions = [Average(self.parset, m, p['avg2'], prefix='facet',
            direction=d, index=2) for d, m in zip(d_list, shifted_data_mapfiles)]
        avg2_data_mapfiles = self.s.run(actions)

        # concatenate all phase-shifted, averaged bands together. Do this twice,
        # once each of the two averaged file sets
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
        for dm, cdm, d in zip(concat_data_mapfiles, concat_corrdata_mapfiles, d_list):
            concat_data_file = read_mapfile(dm)[0]
            concat_corrdata_file = read_mapfile(cdm)[0]
            d.concat_file = concat_data_file
            copy_column(d.concat_file, p['copy']['incol'], p['copy']['outcol'],
                ms_from=concat_corrdata_file)


class FacetSelfcal(Operation):
    """
    Operation to selfcal one or more directions
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetSelfcal, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetSelfcal')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import Average, Concatenate
        from factor.actions.calibrations import Apply, Solve
        from factor.actions.images import MakeImage
        from factor.actions.models import MakeSkymodelFromModelImage, FFT
        from factor.actions.solutions import Smooth, ResetPhases
        from factor.lib.operation_lib import copy_column, make_chunks, merge_chunks
        from factor.operations.hardcoded_param import facet_selfcal as p
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        d_list = self.direction
        if type(d_list) is not list:
            d_list = [d_list]
        bands = self.bands

        # Check state for each direction
        all_done = False
        for i, d in enumerate(d_list):
            if os.path.exists(self.statebasename[i]+'.done'):
                final_parmdb_datamap = os.path.join(self.parset['dir_working'],
                    'datamaps/FacetSelfcal/merged_parmdb_final_{0}.datamap'.
                    format(d.name))
                file = read_mapfile(final_parmdb_datamap)[0]
                d.dirdepparmdb = file
                all_done = True
        if all_done:
            return

        # Make initial data maps for the averaged, phase-shifted datasets
        facet_data_mapfiles = []
        for d in d_list:
            facet_data_mapfiles.append(write_mapfile([d.concat_file], self.name,
            prefix='shifted_vis', working_dir=self.parset['dir_working']))

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
        actions = [MakeImage(self.parset, m, p['imager0'], prefix='facet_selfcal0',
            direction=d) for d, m in zip(d_list, avg_data_mapfiles)]
        image0_basenames_mapfiles = self.s.run(actions)

        if not self.parset['use_ftw']:
            self.log.info('Making sky model (facet model #0)...')
            actions = [MakeSkymodelFromModelImage(self.parset, m, p['model0'],
                prefix='facet_selfcal0', direction=d) for d, m in zip(d_list,
                image0_basenames_mapfiles)]
            skymodels0_mapfiles = self.s.run(actions)

        # Chunk phase-shifted concatenated MS of all bands in time. Chunks should
        # be the length of the desired amplitude solution interval (otherwise they
        # must be concatenated together later before amp. calibration)
        self.log.info('Dividing dataset into chunks...')
        chunks_list = []
        for d in d_list:
            chunks_list.append(make_chunks(d.concat_file, d.solint_a,
            	self.parset, 'facet_chunk', direction=d))
        chunk_data_mapfiles = []
        chunk_parmdb_mapfiles = []
        chunk_model_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_data_mapfiles.append(write_mapfile([chunk.file for chunk in chunks],
                self.name, prefix='chunks_vis', working_dir=self.parset['dir_working'],
                direction=d_list[i]))
            if self.parset['use_ftw']:
                model0 = read_mapfile(image0_basenames_mapfiles[i])
                chunk_model_mapfiles.append(write_mapfile(model0*len(chunks),
                    self.name, prefix='chunks_model0', working_dir=self.parset['dir_working'],
                    direction=d_list[i]))
            else:
                skymodel0 = read_mapfile(skymodels0_mapfiles[i])
                chunk_model_mapfiles.append(write_mapfile(skymodel0*len(chunks),
                    self.name, prefix='chunks_skymodel0', working_dir=self.parset['dir_working'],
                    direction=d_list[i]))
            chunk_parmdb_mapfiles.append(write_mapfile([chunk.parmdb_phaseonly1 for chunk in chunks],
                self.name, prefix='chunk_parmdb_phaseonly1', working_dir=self.parset['dir_working'],
                direction=d_list[i]))

        if self.parset['use_ftw']:
            self.log.debug('FFTing model image (facet model #0)...')
            actions = [FFT(self.parset, dm, mm, p['model0'], prefix='fft0',
            	direction=d, index=0) for d, dm, mm in zip(d_list,
            	chunk_data_mapfiles, chunk_model_mapfiles)]
            self.s.run(actions)

        self.log.info('Solving for phase solutions and applying them (#1)...')
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseonly1'].copy()
            p_d['timestep'] = d.solint_p
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=mm,
            prefix='facet_phaseonly', direction=d, index=0)
            for d, dm, pd, mm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_model_mapfiles)]
        self.s.run(actions)

        self.log.info('Imaging (facet image #1)...')
        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            merged_data_mapfiles.append(write_mapfile([merge_chunks(
                [chunk.file for chunk in chunks], prefix='image1')], self.name, prefix='merged_vis',
                working_dir=self.parset['dir_working'], direction=d_list[i]))

        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg1'], prefix='facet',
            direction=d, index=1) for d, m in zip(d_list, merged_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Imaging...')
        actions = [MakeImage(self.parset, m, p['imager1'], prefix='facet_selfcal1',
            direction=d) for d, m in zip(d_list, avg_data_mapfiles)]
        image1_basenames_mapfiles = self.s.run(actions)

        if not self.parset['use_ftw']:
            self.log.info('Making sky model (facet model #1)...')
            actions = [MakeSkymodelFromModelImage(self.parset, m, p['model1'],
                prefix='facet_selfcal1', direction=d) for d, m in zip(d_list,
                image1_basenames_mapfiles)]
            skymodels1_mapfiles = self.s.run(actions)

        chunk_parmdb_mapfiles = []
        chunk_model_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            if self.parset['use_ftw']:
                model1 = read_mapfile(image1_basenames_mapfiles[i])
                chunk_model_mapfiles.append(write_mapfile(model1*len(chunks),
                    self.name, prefix='chunks_model1', working_dir=self.parset['dir_working'],
                    direction=d_list[i]))
            else:
                skymodel1 = read_mapfile(skymodels1_mapfiles[i])
                chunk_model_mapfiles.append(write_mapfile(skymodel1*len(chunks),
                    self.name, prefix='chunks_skymodel1', working_dir=
                    self.parset['dir_working'], direction=d_list[i]))
            chunk_parmdb_mapfiles.append(write_mapfile([chunk.parmdb_phaseonly2 for chunk in chunks],
                self.name, prefix='chunk_parmdb_phaseonly2', working_dir=
                self.parset['dir_working'], direction=d_list[i]))

        if self.parset['use_ftw']:
            self.log.debug('FFTing model image (facet model #1)...')
            actions = [FFT(self.parset, dm, mm, p['model1'], prefix='fft1',
            	direction=d, index=1) for d, dm, mm in zip(d_list,
            	chunk_data_mapfiles, chunk_model_mapfiles)]
            self.s.run(actions)

        self.log.info('Solving for phase solutions and applying them (#2)...')
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseonly2'].copy()
            p_d['timestep'] = d.solint_p
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=mm,
            parmdb_datamap=pm, prefix='facet_phaseonly', direction=d, index=1)
            for d, dm, pd, mm, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_model_mapfiles, chunk_parmdb_mapfiles)]
        self.s.run(actions)

        self.log.info('Imaging (facet image #2)...')
        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for chunks in chunks_list:
            merged_data_mapfiles.append(write_mapfile([merge_chunks(
                [chunk.file for chunk in chunks], prefix='image2')], self.name, prefix='merged_vis',
                working_dir=self.parset['dir_working']))

        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg2'], prefix='facet',
            direction=d, index=2) for d, m in zip(d_list, merged_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Imaging...')
        actions = [MakeImage(self.parset, m, p['imager2'], prefix='facet_selfcal2',
            direction=d) for d, m in zip(d_list, avg_data_mapfiles)]
        image2_basenames_mapfiles = self.s.run(actions)

        if not self.parset['use_ftw']:
            self.log.info('Making sky model (facet model #2)...')
            actions = [MakeSkymodelFromModelImage(self.parset, m, p['model2'],
                prefix='facet_selfcal2', direction=d) for d, m in zip(d_list,
                image2_basenames_mapfiles)]
            skymodels2_mapfiles = self.s.run(actions)

        chunk_parmdb_phaseamp_phase1_mapfiles = []
        chunk_parmdb_phaseamp_amp1_mapfiles = []
        chunk_model_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            if self.parset['use_ftw']:
                model2 = read_mapfile(image2_basenames_mapfiles[i])
                chunk_model_mapfiles.append(write_mapfile(model2*len(chunks),
                    self.name, prefix='chunks_model2', working_dir=self.parset['dir_working'],
                    direction=d_list[i]))
            else:
                skymodel2 = read_mapfile(skymodels2_mapfiles[i])
                chunk_model_mapfiles.append(write_mapfile(skymodel2*len(chunks),
                    self.name, prefix='chunks_skymodel2', working_dir=self.parset['dir_working']))
            chunk_parmdb_phaseamp_phase1_mapfiles.append(write_mapfile(
                [chunk.parmdb_phaseamp_phase1 for chunk in chunks],
                self.name, prefix='chunk_parmdb_phase1', working_dir=
                self.parset['dir_working'], direction=d_list[i]))
            chunk_parmdb_phaseamp_amp1_mapfiles.append(write_mapfile(
                [chunk.parmdb_phaseamp_amp1 for chunk in chunks],
                self.name, prefix='chunk_parmdb_amp1', working_dir=
                self.parset['dir_working'], direction=d_list[i]))

        if self.parset['use_ftw']:
            self.log.debug('FFTing model image (facet model #2)...')
            actions = [FFT(self.parset, dm, mm, p['model2'], prefix='fft2',
            	direction=d, index=2) for d, dm, mm in zip(d_list,
            	chunk_data_mapfiles, chunk_model_mapfiles)]
            self.s.run(actions)

        self.log.info('Solving for amplitude solutions and applying them (#1)...')
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseamp1_phaseonly'].copy()
            p_d['timestep'] = d.solint_p
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=mm,
            parmdb_datamap=pm, prefix='facet_phaseonly', direction=d, index=2)
            for d, dm, pd, mm, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_model_mapfiles, chunk_parmdb_phaseamp_phase1_mapfiles)]
        self.s.run(actions)
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseamp1_amponly'].copy()
            p_d['timestep'] = d.solint_a
            p_d['chunksize'] = d.solint_a
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=mm,
            parmdb_datamap=pm, prefix='facet_amponly', direction=d, index=2)
            for d, dm, pd, mm, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_model_mapfiles, chunk_parmdb_phaseamp_amp1_mapfiles)]
        self.s.run(actions)

        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for chunks in chunks_list:
            merged_data_mapfiles.append(write_mapfile([merge_chunks(
                [chunk.file for chunk in chunks], prefix='image3')], self.name, prefix='merged_vis',
                working_dir=self.parset['dir_working'], direction=d_list[i]))

        self.log.info('Merging instrument parmdbs...')
        merged_parmdb_phaseamp_amp1_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            concat_file = read_mapfile(merged_data_mapfiles[i])
            merged_parmdb_phaseamp_amp1_mapfiles.append(write_mapfile(
                [self.merge_chunk_parmdbs([chunk.parmdb_phaseamp_amp1 for chunk
                in chunks], concat_file, prefix='merged_amps1')], self.name, prefix='merged_amps1',
                working_dir=self.parset['dir_working']))

        self.log.info('Smoothing amplitude solutions...')
        actions = [Smooth(self.parset, dm, p['smooth_amp1'], pm,
            prefix='facet_amp', direction=d, index=2)
            for d, dm, pm in zip(d_list, merged_data_mapfiles,
            merged_parmdb_phaseamp_amp1_mapfiles)]
        self.s.run(actions)

        self.log.info('Applying amplitude solutions...')
        actions = [Apply(self.parset, dm, p['apply_amp1'],
            pm, prefix='facet_amp', direction=d, index=2) for d, dm, pm in
            zip(d_list, merged_data_mapfiles,
            merged_parmdb_phaseamp_amp1_mapfiles)]
        self.s.run(actions)

        self.log.info('Imaging (facet image #3)...')
        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg3'], prefix='facet',
            direction=d, index=3) for d, m in zip(d_list, merged_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Imaging...')
        actions = [MakeImage(self.parset, m, p['imager3'], prefix='facet_selfcal3',
            direction=d) for d, m in zip(d_list, avg_data_mapfiles)]
        image3_basenames_mapfiles = self.s.run(actions)

        if not self.parset['use_ftw']:
            self.log.info('Making sky model (facet model #3)...')
            actions = [MakeSkymodelFromModelImage(self.parset, m, p['model3'],
                prefix='facet_selfcal3', direction=d) for d, m in zip(d_list,
                image3_basenames_mapfiles)]
            skymodels3_mapfiles = self.s.run(actions)

        self.log.debug('Resetting phases...')
        actions = [ResetPhases(self.parset, dm, p['reset_phases'], pm,
            prefix='facet', direction=d)
            for d, dm, pm in zip(d_list, merged_data_mapfiles,
            merged_parmdb_phaseamp_amp1_mapfiles)]
        self.s.run(actions)

        self.log.info('Preapplying amplitude solutions...')
        chunk_parmdb_phaseamp_amp1_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            chunk_parmdb_phaseamp_amp1_mapfiles.append(write_mapfile(
                [merged_parmdb_phaseamp_amp1_mapfiles[i]]*len(chunks),
                self.name, prefix='chunk_parmdb_amp1_smoothed', working_dir=
                self.parset['dir_working'], direction=d_list[i]))
        actions = [Apply(self.parset, dm, p['apply_amp2'],
            pm, prefix='facet_amp', direction=d, index=2) for d, dm, pm in
            zip(d_list, chunk_data_mapfiles,
            chunk_parmdb_phaseamp_amp1_mapfiles)]
        self.s.run(actions)

        self.log.info('Solving for amplitude solutions (#2)...')
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseamp2_phaseonly'].copy()
            p_d['timestep'] = d.solint_p
            p_list.append(p_d)
        chunk_parmdb_phaseamp_phase2_mapfiles = []
        chunk_parmdb_phaseamp_amp2_mapfiles = []
        chunk_model_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            if self.parset['use_ftw']:
                model3 = read_mapfile(image3_basenames_mapfiles[i])
                chunk_model_mapfiles.append(write_mapfile(model3*len(chunks),
                    self.name, prefix='chunks_model3', working_dir=self.parset['dir_working'],
                    direction=d_list[i]))
            else:
                skymodel3 = read_mapfile(skymodels3_mapfiles[i])
                chunk_model_mapfiles.append(write_mapfile(skymodel3*len(chunks),
                    self.name, prefix='chunks_skymodel3', working_dir=
                    self.parset['dir_working'], direction=d_list[i]))
            chunk_parmdb_phaseamp_phase2_mapfiles.append(write_mapfile(
                [chunk.parmdb_phaseamp_phase2 for chunk in chunks],
                self.name, prefix='chunk_parmdb_phase2', working_dir=
                self.parset['dir_working'], direction=d_list[i]))
            chunk_parmdb_phaseamp_amp2_mapfiles.append(write_mapfile(
                [chunk.parmdb_phaseamp_amp2 for chunk in chunks],
                self.name, prefix='chunk_parmdb_amp2', working_dir=
                self.parset['dir_working'], direction=d_list[i]))

        if self.parset['use_ftw']:
            self.log.debug('FFTing model image (facet model #2)...')
            actions = [FFT(self.parset, dm, mm, p['model3'], prefix='fft3',
            	direction=d, index=3) for d, dm, mm in zip(d_list,
            	chunk_data_mapfiles, chunk_model_mapfiles)]
            self.s.run(actions)

        actions = [Solve(self.parset, dm, pd, model_datamap=mm,
            parmdb_datamap=pm, prefix='facet_phaseonly', direction=d, index=3)
            for d, dm, pd, mm, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_model_mapfiles, chunk_parmdb_phaseamp_phase2_mapfiles)]
        self.s.run(actions)
        p_list = []
        for d in d_list:
            p_d = p['solve_phaseamp2_amponly'].copy()
            p_d['timestep'] = d.solint_a
            p_d['chunksize'] = d.solint_a
            p_list.append(p_d)
        actions = [Solve(self.parset, dm, pd, model_datamap=mm,
            parmdb_datamap=pm, prefix='facet_amponly', direction=d, index=3)
            for d, dm, pd, mm, pm in zip(d_list, chunk_data_mapfiles, p_list,
            chunk_model_mapfiles, chunk_parmdb_phaseamp_amp2_mapfiles)]
        self.s.run(actions)

        self.log.debug('Merging chunks...')
        merged_data_mapfiles = []
        for chunks in chunks_list:
            merged_data_mapfiles.append(write_mapfile([merge_chunks(
            	[chunk.file for chunk in chunks], prefix='image4')],
            	self.name, prefix='merged_vis',
            	working_dir=self.parset['dir_working'], direction=d_list[i]))

        self.log.info('Merging instrument parmdbs...')
        merged_parmdb_phaseamp_amp2_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            concat_file = read_mapfile(merged_data_mapfiles[i])
            merged_parmdb_phaseamp_amp2_mapfiles.append(write_mapfile(
            	[self.merge_chunk_parmdbs([chunk.parmdb_phaseamp_amp2 for
            	chunk in chunks], concat_file, prefix='merged_parmdb_amps2')],
            	self.name, prefix='merged_amps2',
            	working_dir=self.parset['dir_working'], direction=d_list[i]))
        merged_parmdb_phaseamp_phase2_mapfiles = []
        for i, chunks in enumerate(chunks_list):
            concat_file = read_mapfile(merged_data_mapfiles[i])
            merged_parmdb_phaseamp_phase2_mapfiles.append(write_mapfile(
            	[self.merge_chunk_parmdbs([chunk.parmdb_phaseamp_phase2 for
            	chunk in chunks], concat_file,
            	prefix='merged_parmdb_phases2')], self.name,
            	prefix='merged_phases2',
            	working_dir=self.parset['dir_working'], direction=d_list[i]))

        self.log.info('Smoothing amplitude solutions...')
        actions = [Smooth(self.parset, dm, p['smooth_amp2'], pm,
            prefix='facet_amp', direction=d, index=3)
            for d, dm, pm in zip(d_list, merged_data_mapfiles,
            merged_parmdb_phaseamp_amp2_mapfiles)]
        self.s.run(actions)

        self.log.info('Applying amplitude solutions...')
        actions = [Apply(self.parset, dm, p['apply_amp3'],
            pm, prefix='facet_amp', direction=d, index=3) for d, dm, pm in
            zip(d_list, merged_data_mapfiles,
            merged_parmdb_phaseamp_amp2_mapfiles)]
        self.s.run(actions)

        self.log.info('Imaging (facet image #4)...')
        self.log.debug('Averaging in preparation for imaging...')
        actions = [Average(self.parset, m, p['avg4'], prefix='facet',
            direction=d, index=4) for d, m in zip(d_list, merged_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.debug('Imaging...')
        actions = [MakeImage(self.parset, m, p['imager4'], prefix='facet_selfcal4',
            direction=d) for d, m in zip(d_list, avg_data_mapfiles)]
        image3_basenames_mapfiles = self.s.run(actions)

        # Check images if interactive=True
        # If not OK: stop and reset state for this direction
        # If OK, continue
        if self.parset['interactive']:
            import pyrap.images as pim
            for d in d_list:
                im = pim.image([image1[0], image2[0], image3[0], image4[0]])
                im.view()
                prompt = "Continue processing (y/n)? "
                answ = raw_input(prompt)
                while answ.lower() not in  ['y', 'n', 'yes', 'no']:
                    answ = raw_input(prompt)
                if answ.lower() in ['n', 'no']:
                    self.log.info('Resetting facet {0} and exiting...'.format(d.name))
                    sys.exit()

        self.log.info('Merging final instrument parmdbs...')
        merged_parmdb_final_mapfiles = []
        for i, d in enumerate(d_list):
            phases2_final = read_mapfile(merged_parmdb_phaseamp_phase2_mapfiles[i])[0]
            smoothed_amps1_final = read_mapfile(merged_parmdb_phaseamp_amp1_mapfiles[i])[0]
            smoothed_amps2_final = read_mapfile(merged_parmdb_phaseamp_amp2_mapfiles[i])[0]
            concat_file = read_mapfile(merged_data_mapfiles[i])[0]
            merged_parmdb_final_mapfiles.append(write_mapfile(
                [self.merge_parmdbs(phases2_final, smoothed_amps1_final, smoothed_amps2_final,
                concat_file)], self.name, prefix='merged_parmdb_final',
                working_dir=self.parset['dir_working'], direction=d))

        self.log.info('Smoothing amplitude solutions...')
        actions = [Smooth(self.parset, dm, p['smooth_amp3'], pm,
            prefix='facet_amp', direction=d, index=4)
            for d, dm, pm in zip(d_list, merged_data_mapfiles,
            merged_parmdb_final_mapfiles)]
        self.s.run(actions)

        # Save files to the direction objects
        for d, m in zip(d_list, merged_parmdb_final_mapfiles):
            f = read_mapfile(m)[0]
            d.dirdepparmdb = f


    def merge_chunk_parmdbs(self, inparmdbs, concat_file, prefix='merged',
        clobber=False):
        """Merges parmdbs"""
        import os
        import lofar.parmdb
        import pyrap.tables as pt

        root_dir = inparmdbs[0].split('chunks')[0]
        outparmdb = '{0}/{1}_instrument'.format(root_dir, prefix)
        if os.path.exists(outparmdb):
            if clobber:
                os.system('rm -rf {0}'.format(outparmdb))
            else:
                return outparmdb

        os.system('cp -r {0} {1}'.format(inparmdbs[0], outparmdb))
        pdb_concat = lofar.parmdb.parmdb(outparmdb)

        for parmname in pdb_concat.getNames():
            pdb_concat.deleteValues(parmname)

        for inparmdb in inparmdbs:
            pdb = lofar.parmdb.parmdb(inparmdb)
            for parmname in pdb.getNames():
                v = pdb.getValuesGrid(parmname)
                pdb_concat.addValues(v)
        pdb_concat.flush()

        return outparmdb


    def merge_parmdbs(self, parmdb_p, pre_apply_parmdb, parmdb_a, ms, clobber=True):
        """Merges amp+phase parmdbs"""
        import lofar.parmdb
        import pyrap.tables as pt
        import numpy as np

        parmdb_out = parmdb_p.split('phases2')[0] + 'final' + parmdb_p.split('phases2')[1]
        if os.path.exists(parmdb_out):
            if clobber:
                os.system('rm -rf {0}'.format(parmdb_out))
            else:
                return parmdb_out
        pdb_out = lofar.parmdb.parmdb(parmdb_out, create=True)

        pol_list = ['0:0', '1:1']
        gain = 'Gain'
        anttab = pt.table(ms + '::ANTENNA', ack=False)
        antenna_list = anttab.getcol('NAME')
        anttab.close()

        # Copy over the CommonScalar phases and TEC
        pdb_p   = lofar.parmdb.parmdb(parmdb_p)
        for parmname in pdb_p.getNames():
            parms = pdb_p.getValuesGrid(parmname)
            pdb_out.addValues(parms)

        # Get amplitude solutions
        pdb_pre = lofar.parmdb.parmdb(pre_apply_parmdb)
        pdb_a = lofar.parmdb.parmdb(parmdb_a)
        parms_pre = pdb_pre.getValuesGrid("*")
        parms_a = pdb_a.getValuesGrid("*")

        # Get array sizes and initialize values using first antenna (all
        # antennas should be the same)
        parmname = 'Gain:0:0:Real:' + antenna_list[0]
        N_times, N_freqs = parms_a[parmname]['values'].shape
        times = parms_a[parmname]['times'].copy()
        timewidths = parms_a[parmname]['timewidths'].copy()
        freqs = parms_a[parmname]['freqs'].copy()
        freqwidths = parms_a[parmname]['freqwidths'].copy()
        parms = {}
        v = {}
        v['times'] = times
        v['timewidths'] = timewidths
        v['freqs'] = freqs
        v['freqwidths'] = freqwidths

        # Multiply gains
        for pol in pol_list:
            for antenna in antenna_list:
                real1 = np.copy(parms_pre[gain+':'+pol+':Real:'+antenna]['values'])
                real2 = np.copy(parms_a[gain+':' +pol+':Real:'+antenna]['values'])
                imag1 = np.copy(parms_pre[gain+':'+pol+':Imag:'+antenna]['values'])
                imag2 = np.copy(parms_a[gain+':'+pol+':Imag:'+antenna]['values'])

                G1 = real1 + 1j * imag1
                G2 = real2 + 1j * imag2
                Gnew = G1 * G2
                v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)

                parmname = gain + ':' + pol + ':Imag:' + antenna
                v['values'] = np.copy(np.imag(Gnew))
                parms[parmname] = v.copy()

                parmname = gain + ':' + pol + ':Real:' + antenna
                v['values'] = np.copy(np.real(Gnew))
                parms[parmname] = v.copy()

        pdb_out.addValues(parms)
        pdb_out.flush()

        return parmdb_out


class FacetAddAll(Operation):
    """
    Operation to add all sources in the facet to data
    """
    def __init__(self, parset, bands, direction=None, reset=False,
        name='FacetAddAll'):
        super(FacetAddAll, self).__init__(parset, bands, direction=direction,
            reset=reset, name=name)


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import PhaseShift
        from factor.actions.models import MakeFacetSkymodel
        from factor.actions.calibrations import Add
        from factor.operations.hardcoded_param import facet_add_all as p
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        d = self.direction
        bands = self.bands

        # Check state
        if os.path.exists(self.statebasename+'.done'):
            shifted_data_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/FacetAddAll/PhaseShift/facet_output_{0}.datamap'.
                format(d.name))
            files = read_mapfile(shifted_data_mapfile)
            for band, f in zip(bands, files):
                band.shifted_data_file = f
            return

        # Make initial data maps for the empty datasets, their dir-indep
        # instrument parmdbs, and their dir-indep sky models
        subtracted_all_mapfile = write_mapfile([band.file for band in bands],
            self.name, prefix='subtracted_all', working_dir=self.parset['dir_working'])
        dir_indep_parmdbs_mapfile = write_mapfile([band.dirindparmdb for
        	band in bands], self.name, prefix='dir_indep_parmdbs',
        	working_dir=self.parset['dir_working'])
        dir_indep_skymodels_mapfile = write_mapfile([band.skymodel_dirindep
            for band in bands], self.name, prefix='dir_indep_skymodels',
            working_dir=self.parset['dir_working'])

        self.log.info('Selecting sources for this direction...')
        action = MakeFacetSkymodel(self.parset, dir_indep_skymodels_mapfile,
            p['select'], d, prefix='all_final', cal_only=False)
        dir_indep_all_skymodels_mapfile = action.run()

        self.log.info('Adding sources for this direction...')
        self.parset['use_ftw'] = False
        action = Add(self.parset, subtracted_all_mapfile, p['add'],
            dir_indep_all_skymodels_mapfile, dir_indep_parmdbs_mapfile,
            prefix='facet_dirindep', direction=d)
        action.run()

        self.log.info('Phase shifting DATA...')
        action = PhaseShift(self.parset, subtracted_all_mapfile, p['shift'],
            prefix='facet', direction=d)
        shifted_data_mapfile = action.run()

        # Save files to the band objects
        files = read_mapfile(shifted_data_mapfile)
        for band, f in zip(bands, files):
            band.shifted_data_file = f


class FacetImage(Operation):
    """
    Operation to image the full facet
    """
    def __init__(self, parset, bands, direction=None, reset=False,
        name='FacetImage'):
        super(FacetImage, self).__init__(parset, bands, direction=direction,
            reset=reset, name=name)


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.visibilities import Average, Concatenate
        from factor.actions.calibrations import Apply
        from factor.actions.images import MakeImage
        from factor.actions.models import MakeFacetSkymodel, MakeSkymodelFromModelImage
        from factor.lib.operation_lib import copy_column
        from factor.operations.hardcoded_param import facet_image as p
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        d_list = self.direction
        if type(d_list) is not list:
            d_list = [d_list]
        bands = self.bands

        # Check state for each direction
        all_done = False
        for i, d in enumerate(d_list):
            if os.path.exists(self.statebasename[i]+'.done'):
                model_mapfile = os.path.join(self.parset['dir_working'],
                'datamaps/FacetImage/final_model_{0}.datamap'.
                format(d.name))
                file = read_mapfile(model_mapfile)[0]
                d.skymodel_dirdep = file
                all_done = True
        if all_done:
            return

        # Make initial data maps for the phase-shifted datasets and their dir-dep
        # instrument parmdbs
        shifted_data_mapfiles = []
        dir_dep_parmdbs_mapfiles = []
        for d in d_list:
            shifted_data_mapfiles.append(write_mapfile([band.
            	shifted_data_file for band in bands], self.name,
            	prefix='shifted', working_dir=self.parset['dir_working']))
            dir_dep_parmdbs_mapfiles.append(write_mapfile([d.
            	dirdepparmdb]*len(bands), self.name,
            	prefix='dir_indep_parmdbs',
            	working_dir=self.parset['dir_working']))

        self.log.info('Applying direction-dependent calibration...')
        actions = [Apply(self.parset, dm, p['apply_dirdep'],
            pm, prefix='facet_dirdep', direction=d) for d, dm, pm in zip(d_list,
            shifted_data_mapfiles, dir_dep_parmdbs_mapfiles)]
        self.s.run(actions)

        self.log.info('Averaging DATA...')
        actions = [Average(self.parset, m, p['avg'], prefix='facet',
            direction=d, index=1) for d, m in zip(d_list, shifted_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        self.log.info('Concatenating bands...')
        actions = [Concatenate(self.parset, m, p['concat'],
            prefix='facet_bands', direction=d, index=0) for d, m in zip(d_list,
            avg_data_mapfiles)]
        concat_data_mapfiles = self.s.run(actions)

        # Image with facet region as initial maskmask=[[110,110,150,145],'mycleanbox.txt',
#                                     'myimage.mask','myregion.rgn']
        self.log.info('Imaging...')
        actions = [MakeImage(self.parset, m, p['imager'], prefix='facet_image',
            direction=d) for d, m in zip(d_list, concat_data_mapfiles)]
        image_basenames_mapfiles = self.s.run(actions)

        if self.parset['use_ftw']:
            # Save image basenames to the direction objects
            files = read_mapfile(image_basenames_mapfiles)
            for d, f in zip(d_list, files):
                d.skymodel_dirdep = f
        else:
            self.log.info('Making sky model...')
            actions = [MakeSkymodelFromModelImage(self.parset, m, p['model'],
                prefix='facet_image', direction=d) for d, m in zip(d_list,
                image_basenames_mapfiles)]
            skymodels_mapfiles = self.s.run(actions)

            self.log.info('Selecting only those sources inside the facet...')
            actions = [MakeFacetSkymodel(self.parset, m,
                p['select'], d, prefix='all_final', cal_only=False) for d, m in
                zip(d_list, skymodels_mapfiles)]
            dir_dep_all_skymodels_mapfiles = self.s.run(actions)

            # Save files to the direction objects
            for d, mf in zip(d_list, dir_dep_all_skymodels_mapfiles):
                d.skymodel_dirdep = read_mapfile(mf)[0]

        # Make final data maps
        for d in d_list:
            write_mapfile([d.skymodel_dirdep], self.name,
            	prefix='final_model', working_dir=self.parset['dir_working'])



class FacetSubAll(Operation):
    """
    Operation to subtract final facet sky model from visibilties
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetSubAll, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetSubAll')


    def run_steps(self):
        """
        Run the steps for this operation
        """
        from factor.actions.calibrations import Subtract
        from factor.actions.models import FFT
        from factor.operations.hardcoded_param import facet_sub_all as p
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        d = self.direction
        bands = self.bands

        # Check state
#         if os.path.exists(self.statebasename+'.done'):
#             self.direction.skymodel_dirdep = 'visdata/allbands_concat_{0}.ms'.format(
#                 self.direction.name)
#             return

        # Make initial data maps for the empty datasets, their dir-dep
        # instrument parmdbs, and their dir-dep sky models
        subtracted_all_mapfile = write_mapfile([band.file for band in bands],
        	self.name, prefix='subtracted_all',
        	working_dir=self.parset['dir_working'])
        dir_dep_parmdbs_mapfile = write_mapfile([d.dirdepparmdb]*len(bands),
            self.name, prefix='dir_dep_parmdbs',
        	working_dir=self.parset['dir_working'])
        dir_dep_models_mapfile = write_mapfile([d.skymodel_dirdep]*len(bands),
            self.name, prefix='dir_dep_skymodels',
        	working_dir=self.parset['dir_working'])

        if self.parset['use_ftw']:
            self.log.debug('FFTing model image (facet model final)...')
            action = FFT(self.parset, subtracted_all_mapfile,
            	dir_dep_models_mapfile, p['fft'], prefix='fft', direction=d)
            action.run()

        self.log.info('Subtracting sources for this direction using final model...')
        action = Subtract(self.parset, subtracted_all_mapfile, p['subtract'], dir_dep_models_mapfile, dir_dep_parmdbs_mapfile,
            prefix='facet_dirdep', direction=d)
        action.run()

        # Clean up files for this direction:
        #    - unneeded MS files
        #    - selfcal images
        self.log.info('Cleaning up files for this direction...')
#         os.system('rm -rf visdata/*{0}*'.format(d.name))


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
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        d = self.direction
        bands = self.bands

        # Check state for each direction
#         if os.path.exists(self.statebasename+'.done'):
#             self.direction.concat_file = 'visdata/allbands_concat_{0}.ms'.format(
#                 self.direction.name)
#             return

        # Make initial data maps for the empty datasets, their dir-indep
        # instrument parmdbs, and their dir-indep sky models
        subtracted_all_mapfile = write_mapfile([band.file for band in bands],
        	self.name, prefix='subtracted_all',
        	working_dir=self.parset['dir_working'])
        dir_dep_parmdbs_mapfile = write_mapfile([d.dirdepparmdb]*len(bands),
            self.name, prefix='dir_dep_parmdbs',
        	working_dir=self.parset['dir_working'])
        dir_dep_models_mapfile = write_mapfile([d.skymodel_dirdep]*len(bands),
            self.name, prefix='dir_dep_skymodels',
        	working_dir=self.parset['dir_working'])

        if self.parset['use_ftw']:
            self.log.debug('FFTing model image (facet model final)...')
            action = FFT(self.parset, subtracted_all_mapfile,
            	dir_dep_models_mapfile, p['fft'], prefix='fft', direction=d)
            action.run()

        self.log.info('Adding sources for this direction...')
        action = Add(self.parset, subtracted_all_mapfile, p['add'],
            dir_indep_cal_skymodels_mapfile, dir_indep_parmdbs_mapfile,
            prefix='facet_dirdep', direction=d)
        action.run()

        self.log.info('Phase shifting DATA...')
        action = PhaseShift(self.parset, subtracted_all_mapfile, p['shift'],
            prefix='facet', direction=d)
        shifted_data_mapfile = action.run()

        # Save files to the band objects
        files = read_mapfile(shifted_data_mapfile)
        for band, f in zip(bands, files):
            band.shifted_data_file = f


def FacetImageFinal(FacetImage):
    """
    Operation to make final facet image
    """
    def __init__(self, parset, bands, direction=None, reset=False):
        super(FacetImageFinal, self).__init__(parset, bands, direction=direction,
            reset=reset, name='FacetImageFinal')

