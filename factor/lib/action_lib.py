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


def make_image_basename(input_datamap, direction=None, prefix=None):
    """
    Define a standard name pattern for imaging files

    Parameters
    ----------
    prefix : str, optional
        String to prepend to the image basename. If None, 'image' is used

    """
    from factor.lib.datamap_lib import read_mapfile
    import re
    import os

    if prefix is None:
        prefix = 'image'
        logging.warn('Prefix of imagebasename not selected, using "{0}".'.format(prefix))

    msfiles = read_mapfile(input_datamap)
    image_basenames = []

    for msfile in msfiles:
        msbase = os.path.basename(msfile)

        if direction is not None:
            try:
                dirtxt = direction.name
            except:
                dirtxt = direction
            image_basenames.append('%s-%s_%s' % (prefix, re.sub(r'.MS|.ms', '', msbase), dirtxt))
        else:
            image_basenames.append('%s-%s' % (prefix, re.sub(r'.MS|.ms', '', msbase)))

    return image_basenames
