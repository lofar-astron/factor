"""
Module that holds all model-related actions

Classes
-------
MakeSkymodelFromModelImage : Action
    Makes a sky model from a CASA model image
MakeFacetSkymodel : Action
    Makes a sky model for a facet
MergeSkymodels : Action
    Merges two sky models
FFT : Action
    FFTs a model into an MS

"""

from factor.lib.action import Action
from factor.lib.action_lib import make_parset_basename, make_pipeline_dirname
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class MakeSkymodelFromModelImage(Action):
    """
    Action to make a sky model from a CASA model image
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
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
        super(MakeSkymodelFromModelImage, self).__init__(op_parset,
            'MakeSkymodelFromModelImage')

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


class MakeFacetSkymodel(Action):
    """
    Action to make a sky model for a facet
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
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
        super(MakeFacetSkymodel, self).__init__(op_parset, 'MakeFacetSkymodel')

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if prefix is None:
            prefix = 'make_facet_skymodel'
        self.prefix = prefix
        self.direction = direction
        self.localdir = localdir
        self.clean = clean
        self.index = index


class MergeSkymodels(Action):
    """
    Action to merge two sky models
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
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
        super(MergeSkymodels, self).__init__(op_parset, 'MergeSkymodels')

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if prefix is None:
            prefix = 'make_image'
        self.prefix = prefix
        self.direction = direction
        self.localdir = localdir
        self.clean = clean
        self.index = index


class FFT(Action):
    """
    Action to FFT a model into an MS
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
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
        super(FFT, self).__init__(op_parset, 'FFT')

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if prefix is None:
            prefix = 'make_image'
        self.prefix = prefix
        self.direction = direction
        self.localdir = localdir
        self.clean = clean
        self.index = index
