"""
Some functions called by multiple actions
"""

def makeimagename(ms, prefix, direction = None):
    """
    define a standard name patter for imaging
    """
    import re
    if direction != None:
        return 'img/%s-%s_%s' % (prefix, re.sub(r'.MS|.ms','',ms), direction.name)
    else:
        return 'img/%s-%s' % (prefix, re.sub(r'.MS|.ms','',ms))
