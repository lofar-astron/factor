"""
Module that stores the version and changelog
"""

# Version number
__version__ = '1.3'

# Change log
def changelog():
    """
    FACTOR Changelog.
    ----------------------------------------------------------------------------

    18/07/2017

        - Options have been added to allow the scales used by WSClean during
        clean to be specified (selfcal_multiscale_scales_pixel and
        facet_multiscale_scales_pixel under [imaging])

        - The atrous_do column in the directions file has been renamed to
        mscale_selfcal_do, as this name better describes its purpose

        - Improved handling of flagged solutions during smoothing of the
        amplitude solutions

    14/06/2017 - Version 1.3

        - Updated to use WSClean v2.4 and LSMTool v1.2.0

    25/04/2017

        - Preaveraging is now done in frequency as well as in time. This
        preaveraging generally improves the S/N of the CS slow-gain solutions
        for fainter sources by a factor of ~ 2
        - An option (min_fraction_per_band) has been added that sets the
        minimum allowed unflagged fraction per band

    11/04/2017

        - WSClean's automasking feature is now used during imaging. The old
        image-mask-image sequence is no longer used during self calibration, but
        can still be used during the final, full-bandwidth facet imaging if
        automask_facet_image = False under the [imaging] section of the
        parset. The option update_selfcal_clean_regions has been removed, as it
        no longer applies (PyBDSF masking has been removed from selfcal)
        - Handling of pipeline failure/interruption has been improved

    07/04/2017

        - The combination of flagging ranges specified by the flag_abstime,
        flag_baseline, and flag_freqrange options can now be set with the
        flag_expr option

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
