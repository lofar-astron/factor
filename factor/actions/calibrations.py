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
            prefix = 'run_bbs'
        self.prefix = prefix
        self.direction = direction
        self.localdir = localdir
        self.clean = clean
        self.index = index

        self.ms = ms
        self.statebasename += makestatebasename(self.ms, prefix, direction, index)
        self.logbasename += makestatebasename(self.ms, prefix, direction, index)
        self.parsetbasename += makestatebasename(self.ms, prefix, direction, index)
        self.parsetname = self.parsetbasename + '.parset'
        self.clean = clean
        self.p = p.copy()
        if '/' in parmdb:
            # path is not local to input ms, so copy it to ms
            destfile = '/'.join([self.ms, os.path.basename(parmdb)])
            if os.path.exists(destfile):
                os.system('rm -rf {0}'.format(destfile))
            shutil.copytree(parmdb, destfile)
            self.parmdb = os.path.basename(parmdb)
        else:
            self.parmdb = parmdb
        self.templatename = '{0}_{1}.parset.tpl'.format(prefix, self.name)
        if skymodel is None:
            self.skymodel = self.modelbasename + 'empty.skymodel'
            os.system('touch {0}'.format(self.skymodel))
        else:
            self.skymodel = skymodel

    def get_command(self):
        template_bbs = env.get_template(self.templatename)
        tmp = template_bbs.render(self.p)
        with open(self.parsetname, 'wb') as f:
            f.write(tmp)

        if 'flags' in self.p:
            flags = self.p['flags']
        else:
            flags = ''
        if self.parmdb is not None:
            # Use of --parmdb-name means parmdb path must be local to MS
            self.cmd = 'calibrate-stand-alone {0} --parmdb-name {1} {2} {3} {4}'.format(
                flags, self.parmdb, self.ms, self.parsetname, self.skymodel)
        else:
            self.cmd = 'calibrate-stand-alone {0} {1} {2} {3}'.format(
                flags, self.ms, self.parsetname, self.skymodel)



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
            prefix = 'run_dppp'
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

        # Deal with empty sky models: (Note: if a facet sky model has no sources, we need
        # simply to copy the visibilities). Set the skip flag in the data map, then
        # deal with them as follows:
#         for band in bands:
#              if band.calmodel_dirindep is None:
#                 copy_column(band.file, p['add']['incol'], p['add']['outcol'])


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
