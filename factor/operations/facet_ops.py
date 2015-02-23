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

        if os.path.exists(self.statebasename+'.done'):
            return

        d = self.direction
        bands = self.bands

        # Make initial data maps for the empty datasets, their dir-indep
        # instrument parmdbs, and their dir-indep sky models
        subtracted_all_mapfile = self.make_datamap([band.file for band in bands],
            'initial_subtracted')
        dir_indep_parmdbs_mapfile = self.make_datamap([band.dirindparmdb for band
            in bands], 'dirindep_parmdbs')
        dir_indep_skymodels_mapfile = self.make_datamap([band.skymodel_dirindep
            for band in bands], 'dirindep_skymodels')

        # Add calibrators from the dir-indep sky model for this direction to the
        # visibilities
        self.log.info('Selecting sources for this direction...')
        action = MakeFacetSkymodel(self.name, dir_indep_skymodels_mapfile, d,
            prefix='cal', cal_only=True)
        dir_indep_cal_skymodels_mapfile = action.get_results()

        self.log.info('Adding sources for this direction...')
        action = Add(self.name, [subtracted_all_mapfile,
            dir_indep_cal_skymodels_mapfile, dir_indep_parmdbs_mapfile], p['add'],
            prefix='facet_dirindep', direction=d)

        # Phase shift to facet center
        self.log.info('Phase shifting DATA...')
        action = PhaseShift(self.name, subtracted_all_mapfile, p['shift1'],
            prefix='facet', direction=d, index=1)
        shifted_data_mapfile = action.get_results()

        # Save files to band objects
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
        from factor.lib.operation_lib import copy_column
        from factor.operations.hardcoded_param import facet_setup as p

        if os.path.exists(self.statebasename+'.done'):
            self.direction.concat_file = 'visdata/allbands_concat_{0}.ms'.format(
                self.direction.name)
            return

        d = self.direction
        bands = self.bands

        # Make initial data maps for the phase-shifted datasets and their dir-indep
        # instrument parmdbs
        shifted_data_mapfile = self.make_datamap([band.shifted_data_file for band
            in bands], 'facet_setup')
        dir_indep_parmdbs_mapfile = self.make_datamap([band.dirindparmdb for band
            in bands], 'dirindep_parmdbs')

        # apply direction-independent calibration
        self.log.info('Applying direction-independent calibration...')
        action = Apply(self.name, [shifted_data_mapfile, None, dir_indep_parmdbs_mapfile],
            p['apply'], prefix='facet_dirindep', direction=d)

        # average to 1 channel per band. Do this twice, once for DATA and once
        # for CORRECTED_DATA
        self.log.info('Averaging DATA...')
        action = Average(self.name, shifted_data_mapfile, p['avg1'], prefix='facet',
            direction=d, index=1)
        avg_data_mapfile = action.get_results()
        self.log.info('Averaging CORRECTED_DATA...')
        action = Average(self.name, shifted_data_mapfile, p['avg2'], prefix='facet',
            direction=d, index=2)
        avg_corrdata_mapfile = action.get_results()

        # concatenate all phase-shifted, averaged bands together. Do this twice,
        # once each of the two averaged file sets
        self.log.info('Concatenating bands...')
        action = Concatenate(self.name, avg_data_mapfile, p['concat1'],
            prefix='facet_bands', direction=d)
        concat_data_mapfile = action.get_results()
        action = Concatenate(self.name, avg_corrdata_mapfile, p['concat2'],
            prefix='facet_bands', direction=d)
        concat_corrdata_mapfile = action.get_results()

        # Copy over CORRECTED_DATA column from second concat file, so that
        # first concatenated file has both DATA and dir-indep CORRECTED_DATA
        concat_data_file = read_mapfile(concat_data_mapfile)
        concat_corrdata_file = read_mapfile(concat_corrdata_mapfile)
        d.concat_file = 'visdata/allbands_concat_{0}.ms'.format(d.name)
        os.system('cp -r {0} {1}'.format(concat_data_file, d.concat_file))
        copy_column(d.concat_file, p['copy']['incol'], p['copy']['outcol'],
            ms_from=concat_corrdata_file)


