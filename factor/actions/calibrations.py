"""
Module that holds all calibration-related actions

Classes
-------
BBS : Action
    Runs BBS calibrate-stand-alone
DPPP : Action
    Runs DPPP Gaincal
Add : Action
    Adds sources
Apply : Action
    Applies solutions
Solve : Action
    Solves for solutions
Subtract : Action
    Subtracts sources

"""

from factor.lib.action import Action
from factor.lib.action_lib import make_basename
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class BBS(Action):
    """
    Action to run BBS
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None, name='BBS'):
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
        super(BBS, self).__init__(op_parset, name=name)

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


class DPPP(Action):
    """
    Action to run DPPP Gaincal
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


class Add(BBS):
    """
    Action to add sources
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        localdir=None, clean=True, index=None):
        super(Add, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Add')


class Apply(BBS):
    """
    Action to apply solutions
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        localdir=None, clean=True, index=None):
        super(Apply, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Apply')


class Solve(BBS):
    """
    Action to solve for solutions
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        localdir=None, clean=True, index=None):
        super(Solve, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Solve')


class Subtract(BBS):
    """
    Action to subtract sources
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        localdir=None, clean=True, index=None):
        super(Subtract, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Subtract')
