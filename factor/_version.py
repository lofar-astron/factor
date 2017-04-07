"""
Module that stores the version and changelog
"""

# Version number
__version__ = '1.2'

# Change log
def changelog():
    """
    FACTOR Changelog.
    ----------------------------------------------------------------------------

    07/04/2017 - Version 1.2

    28/03/2017

        - An unarchiving tool (unarchivefactor) has been added that can
        unarchive an archive made with archivefactor

    20/03/2017

        - An option (update_selfcal_clean_regions) has been added that
        controls whether user-supplied clean masks are updated during selfcal
        - Intersections in user-supplied clean masks are now detected and an
        error raised

    09/03/2017

        - An archiving tool (archivefactor) has been added that can archive
        the subtracted datasets, the sky models, the instrument tables, the
        selfcal plots, and the calibrated data for one or more directions

    08/03/2017

        - Polynomial sky models generated directly by WSClean during imaging are
        used for prediction, resulting in improved and faster subtraction of
        extended sources.

    24/02/2017 - Version 1.1

        - The subtraction of sources after self calibration has been improved
        - The baseline-dependent averaging used during imaging has been reduced,
        as the previous averaging caused significant smearing away from the
        facet center
        - Many bug fixes

    04/11/2016 - Version 1.0

        - Initial release
    """
    pass
