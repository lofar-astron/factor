"""
Module that holds all solution-related actions

Classes
-------
LoSoTo : Action
    Runs LoSoTo
Smooth : Action
    Smooths and normalizes solutions
ResetPhases : Action
    Resets phases to zero

"""

from factor.lib.action import Action
from factor.lib.action_lib import make_basename
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class LoSoTo(Action):
    """
    Action to run LoSoTo
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None, name='LoSoTo'):
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
        super(LoSoTo, self).__init__(op_parset, name=name)

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


class Smooth(LoSoTo):
    """
    Action to smooth and normalize solutions
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        super(Smooth, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Smooth')


class ResetPhases(LoSoTo):
    """
    Action to reset phases to zero
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        super(ResetPhases, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='ResetPhases')
