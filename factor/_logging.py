#!/usr/bin/env python
# encoding: utf-8
# this modeule sets some propreties of the logging system

import logging

def add_coloring_to_emit_ansi(fn):
    """
    colorize the logging output
    """
    # add methods we need to the class
    def new(*args):
        levelno = args[0].levelno
        if(levelno>=50):
            color = '\x1b[31m' # red
        elif(levelno>=40):
            color = '\x1b[31m' # red
        elif(levelno>=30):
            color = '\x1b[33m' # yellow
        elif(levelno>=20):
            color = '\x1b[32m' # green
        elif(levelno>=10):
            color = '\x1b[35m' # pink
        else:
            color = '\x1b[0m' # normal
        args[0].msg = color + args[0].msg +  '\x1b[0m'  # normal
        return fn(*args)
    return new

def set_level(level):
    """
    Change verbosity
    """
    if level == 'warning':
        logging.root.setLevel(logging.WARNING)
    elif level == 'info':
        logging.root.setLevel(logging.INFO)
    elif level == 'debug':
        logging.root.setLevel(logging.DEBUG)

# logging to file (B/W)
fh = logging.FileHandler('factor.log')
fh.setLevel(logging.DEBUG)
logging.root.addHandler(fh)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)

# logging to console (color)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s - %(message)s')
ch.setFormatter(formatter)
ch.emit =  add_coloring_to_emit_ansi(ch.emit)
logging.root.addHandler(ch)


