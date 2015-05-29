"""
Module that sets up Factor logging
"""
import logging


def add_coloring_to_emit_ansi(fn):
    """
    Colorize the logging output
    """
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
    Change verbosity of console output
    """
    for handler in logging.root.handlers:
        if handler.name == 'console':
            if level == 'warning':
                handler.setLevel(logging.WARNING)
            elif level == 'info':
                handler.setLevel(logging.INFO)
            elif level == 'debug':
                handler.setLevel(logging.DEBUG)


class Whitelist(logging.Filter):
    """
    Filter out any non-Factor loggers
    """
    def filter(self, record):
        if 'factor' in record.name and 'executable_' not in record.name:
            return True
        else:
            return False


def set_log_file(log_file):
    """
    Define and add a file handler
    """
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG) # file always logs everything
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
    fh.setFormatter(formatter)
    fh.addFilter(Whitelist())
    logging.root.addHandler(fh)


# Define and add console handler (in color)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO) # default log level
formatter = logging.Formatter('%(levelname)s - %(name)s - %(message)s')
ch.setFormatter(formatter)
ch.emit = add_coloring_to_emit_ansi(ch.emit)
ch.set_name('console')
ch.addFilter(Whitelist())
logging.root.addHandler(ch)

# Set root level (the handlers will set their own levels)
logging.root.setLevel(logging.DEBUG)
