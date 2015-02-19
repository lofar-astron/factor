"""
Some functions called by multiple actions
"""


def make_basename(prefix, direction=None, index=None):
    """
    Returns a standard name pattern

    Parameters
    ----------
    prefix : str
        A prefix for the name, usually related to the particular operation step
    direction : Direction object or str, optional
        A direction name
    index : int, optional
        An index for the particular operation step
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
