"""
Some functions called by multiple actions
"""
import logging


def make_parset_basename(prefix, direction=None, index=None):
    """
    define a standard name pattern for pipeline parset files
    """
    import re
    import os

    if direction is not None:
        try:
            dirtxt = '_{0}'.format(direction.name)
        except:
            dirtxt = '_{0}'.format(direction)
    else:
        dirtxt = ''
    if index is not None:
        indtxt = '-{0}'.format(index)
    else:
        indtxt = ''

    return '{0}-{1}{2}'.format(prefix, dirtxt, indtxt)


def make_pipeline_dirname(prefix, direction=None, index=None):
    """
    define a standard name pattern for pipeline runtime dir
    """
    import re
    import os

    if direction is not None:
        try:
            dirtxt = '_{0}'.format(direction.name)
        except:
            dirtxt = '_{0}'.format(direction)
    else:
        dirtxt = ''
    if index is not None:
        indtxt = '-{0}'.format(index)
    else:
        indtxt = ''

    return '{0}{1}{2}'.format(prefix, dirtxt, indtxt)
