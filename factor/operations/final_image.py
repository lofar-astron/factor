"""
Operation: final_image
Implements final imaging.
"""

from factor.lib.operation import *
import logging

class final_image(operation):

    def __init__(self, parset):
        operation(parset)
        logging.info('Operation final_image started!')
