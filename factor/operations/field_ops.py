"""
Module that holds all field (non-facet-specific) operations

Classes
-------
FieldMosaic : Operation
    Makes a mosaic of the field from the facet images

"""
import os
import logging
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap

log = logging.getLogger('factor:field_ops')


class FieldMosaic(Operation):
    """
    Operation to mosiac facet images
    """
    def __init__(self, parset, bands, direction, cellsize_arcsec, robust,
        taper_arcsec, min_uv_lambda):
        # Set name from the facet imaging parameters. If the parameters are all
        # the same as those used for selfcal, just use 'FieldMosaic'; otherwise,
        # append the parameters
        if (cellsize_arcsec != parset['imaging_specific']['selfcal_cellsize_arcsec'] or
            robust != parset['imaging_specific']['selfcal_robust'] or
            taper_arcsec != 0.0 or
            min_uv_lambda != parset['imaging_specific']['selfcal_min_uv_lambda']):
            name = 'FieldMosaic_c{0}r{1}t{2}u{3}'.format(round(cellsize_arcsec, 1),
                    round(robust, 2), round(taper_arcsec, 1), round(min_uv_lambda, 1))
        else:
            name = 'FieldMosaic'
        super(MakeMosaic, self).__init__(parset, bands, direction,
            name=name)

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'fieldmosaic_pipeline.parset'.format(infix)

        # Define extra parameters needed for this operation
        input_files = [b.files for b in self.bands]
        input_files_single = []
        for bandfiles in input_files:
            for filename in bandfiles:
                input_files_single.append(filename)
        self.parms_dict.update({'ms_files_single': input_files_single,
                                'ms_files_grouped' : str(input_files) })

    def finalize(self):
        """
        Finalize this operation
        """
        # Delete averaged data as they're no longer needed
        self.direction.cleanup_mapfiles = [
            os.path.join(self.pipeline_mapfile_dir, 'sorted_groups.mapfile_groups')
            ]
        self.log.debug('Cleaning up files (direction: {})'.format(self.direction.name))
        self.direction.cleanup()
