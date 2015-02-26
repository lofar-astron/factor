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
            self.name, prefix='subtracted_all', working_dir=self.parset['dir_working'])
        dir_indep_parmdbs_mapfile = write_mapfile([os.path.join(band.file,
            band.dirindparmdb) for band in bands], self.name,
            prefix='dir_indep_parmdbs', working_dir=self.parset['dir_working'])
        dir_indep_skymodels_mapfile = write_mapfile([band.skymodel_dirindep
            for band in bands], self.name, prefix='dir_indep_parmdbs',
            working_dir=self.parset['dir_working'])

        # Add calibrators from the dir-indep sky model for this direction to the
        # visibilities
        self.log.info('Selecting sources for this direction...')
        action = MakeFacetSkymodel(self.parset, dir_indep_skymodels_mapfile,
            p['model'], d, prefix='cal', cal_only=True)
        dir_indep_cal_skymodels_mapfile = action.run()

        self.log.info('Adding sources for this direction...')
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
    Set up data for selfcal
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

        # Check state for each direction
#         if os.path.exists(self.statebasename+'.done'):
#             self.direction.concat_file = 'visdata/allbands_concat_{0}.ms'.format(
#                 self.direction.name)
#             return

        d_list = self.direction
        if type(d_list) is not list:
            d_list = [d_list]
        bands = self.bands

        # Make initial data maps for the phase-shifted datasets and their dir-indep
        # instrument parmdbs
        shifted_data_mapfiles = []
        dir_indep_parmdbs_mapfiles = []
        for d in d_list:
            shifted_data_mapfiles.append(write_mapfile([band.shifted_data_file for band
                in bands], self.name, prefix='shifted', working_dir=self.parset['dir_working']))
            dir_indep_parmdbs_mapfiles.append(write_mapfile([os.path.join(band.file,
                band.dirindparmdb) for band in bands], self.name,
                prefix='dir_indep_parmdbs', working_dir=self.parset['dir_working']))

        # average to 1 channel per band. Do this twice, once for DATA and once
        # for CORRECTED_DATA
        self.log.info('Averaging DATA...')
        actions = [Average(self.parset, m, p['avg1'], prefix='facet',
            direction=d, index=1) for d, m in zip(d_list, shifted_data_mapfiles)]
        avg_data_mapfiles = self.s.run(actions)

        # apply direction-independent calibration
        self.log.info('Applying direction-independent calibration...')
        actions = [Apply(self.parset, dm, p['apply'],
            pm, prefix='facet_dirindep', direction=d) for d, dm, pm in zip(d_list,
            shifted_data_mapfile, dir_indep_parmdbs_mapfile)]
        self.s.run(actions)

        # concatenate all phase-shifted, averaged bands together. Do this twice,
        # once each of the two averaged file sets
        self.log.info('Concatenating bands...')
        actions = [Concatenate(self.name, m, p['concat1'],
            prefix='facet_bands', direction=d) for d, m in zip(d_list,
            avg_data_mapfile)]
        concat_data_mapfiles = self.s.run(actions)
        actions = [Concatenate(self.name, m, p['concat2'],
            prefix='facet_bands', direction=d) for d, m in zip(d_list,
            avg_data_mapfile)]
        concat_corrdata_mapfiles = self.s.run(actions)

        # Copy over CORRECTED_DATA column from second concat file, so that
        # first concatenated file has both DATA and dir-indep CORRECTED_DATA
        for dm, cdm, d in zip(concat_data_mapfiles, concat_corrdata_mapfile, d_list):
            concat_data_file = read_mapfile(dm)
            concat_corrdata_file = read_mapfile(cdm)
            d.concat_file = concat_data_file
            copy_column(d.concat_file, p['copy']['incol'], p['copy']['outcol'],
                ms_from=concat_corrdata_file)


