"""
Module that holds all visibility-related actions

Classes
-------
DPPP : Action
    Runs DPPP
Avg : Action
    Averages visibilities
PhaseShift : Action
    Phase shifts visibilities
Concat : Action
    Concatenates visibilities

"""

from factor.lib.action import Action
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class DPPP(Action):
    """
    Action to run DPPP
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None, name='DPPP'):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map for CASA model image
        p : dict
            Input parset dict defining model and pipeline parameters
        prefix : str, optional
            Prefix to use for model names
        direction : Direction object, optional
            Direction for this model
        clean : bool, optional
            Remove unneeded files?
        index : int, optional
            Index of action

        """
        super(DPPP, self).__init__(op_parset, name=name)

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if prefix is None:
            prefix = 'make_skymodel'
        self.prefix = prefix
        self.direction = direction
        self.localdir = localdir
        self.clean = clean
        self.index = index


class Avg(DPPP):
    """
    Action to average visibilities
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        localdir=None, clean=True, index=None):
        super(Avg, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Avg')


class PhaseShift(DPPP):
    """
    Action to phase shift visibilities
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        localdir=None, clean=True, index=None):
        super(PhaseShift, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='PhaseShift')


class Concat(DPPP):
    """
    Action to concatenate visibilities
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        localdir=None, clean=True, index=None):
        super(Concat, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Concat')

