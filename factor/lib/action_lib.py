"""
Some functions called by multiple actions
"""
import logging

def makeimagebasename(ms, prefix = None, direction = None):
    """
    define a standard name pattern for imaging
    """
    import re
    if prefix == None: 
        logging.error('Prefix ofr imagebasename not selected, using xxx.')
        prefix = 'xxx'
        
    if direction != None:
        return 'img/%s-%s_%s' % (prefix, re.sub(r'.MS|.ms','',ms), direction.name)
    else:
        return 'img/%s-%s' % (prefix, re.sub(r'.MS|.ms','',ms))
