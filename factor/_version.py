"""
Module that stores the version and changelog
"""

# Version number
__version__ = '1.1'

# Change log
def changelog():
    """
    FACTOR Changelog.
    ----------------------------------------------------------------------------

    08/03/2017

        - Use polynomial sky models generated directly by WSClean during
        imaging, resulting in improved and faster subtraction of extended
        sources.

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
