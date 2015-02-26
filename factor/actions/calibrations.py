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
    def __init__(self, op_parset, vis_datamap, p, model_datamap=None,
        parmdb_datamap=None, prefix=None, direction=None, clean=True,
        index=None, name='BBS'):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        vis_datamap : data map
            Input data map for MS files
        p : dict
            Input parset dict defining model and pipeline parameters
        model_datamap : data map
            Input data map for sky model files
        parmdb_datamap : data map
            Input data map for parmdb files
        prefix : str, optional
            Prefix to use for run
        direction : Direction object, optional
            Direction for this action
        clean : bool, optional
            Remove unneeded files?
        index : int, optional
            Index of action

        """
        super(BBS, self).__init__(op_parset, name, prefix=prefix,
            direction=direction, index=index)

        # Store input parameters
        self.vis_datamap = vis_datamap
        self.model_datamap = model_datamap
        self.parmdb_datamap = parmdb_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'run_bbs'
        self.clean = clean
        factor_working_dir = op_parset['dir_working']
        if self.direction is None:
            self.working_dir = '{0}/visdata/{1}/'.format(factor_working_dir,
                self.op_name)
        else:
            self.working_dir = '{0}/visdata/{1}/{2}/'.format(factor_working_dir,
                self.op_name, self.direction.name)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        # Define parset name
        self.parset_file = self.parsetbasename + '.parset'
        self.templatename = '{0}_{1}.parset.tpl'.format(prefix, self.name.lower())


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import write_mapfile

        self.p['vis_datamap'] = self.vis_datamap
        if self.model_datamap is not None:
            self.p['skymodel_datamap'] = self.model_datamap
        self.p['parmdb_datamap'] = self.parmdb_datamap


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        self.p['parset'] = os.path.abspath(self.parset_file)
        self.p['outputdir'] = os.path.abspath(self.working_dir)
        self.p['lofarroot'] = self.op_parset['lofarroot']
        self.p['parset'] = self.parset_file
        if 'flags' not in self.p:
            self.p['flags'] = ''

        if self.model_datamap is not None:
            template = env.get_template('bbs.pipeline.parset.tpl')
        else:
            template = env.get_template('bbs_nomodel.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        template = env.get_template(self.templatename)
        tmp = template.render(self.p)
        with open(self.parset_file, 'w') as f:
            f.write(tmp)


    def make_pipeline_config_parset(self):
        """
        Writes the pipeline configuration parset
        """
        template = env.get_template('pipeline.cfg.tpl')
        tmp = template.render(self.op_parset)
        with open(self.pipeline_config_file, 'w') as f:
            f.write(tmp)


    def get_results(self):
        """
        Return results
        """
        return None


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
    def __init__(self, op_parset, vis_datamap, p, model_datamap,
        parmdb_datamap, prefix=None, direction=None, clean=True,
        index=None):
        super(Add, self).__init__(op_parset, vis_datamap, p,
            model_datamap=model_datamap, parmdb_datamap=parmdb_datamap,
            prefix=prefix, direction=direction, clean=clean, index=index,
            name='Add')

        # Deal with empty sky models: (Note: if a facet sky model has no sources, we need
        # simply to copy the visibilities).
        from factor.lib.datamap_lib import read_mapfile, read_mapfile_flags, set_mapfile_flags
        from factor.lib.operation_lib import copy_column

        model_files = read_mapfile(self.model_datamap)
        model_flags = read_mapfile_flags(self.model_datamap)
        vis_files = read_mapfile(self.vis_datamap)
        for model_file, model_flag, vis_file in zip(model_files, model_flags, vis_files):
            if model_flag:
                self.log.info('Skipping add for {0}'.format(vis_file))
                copy_column(vis_file, self.p['incol'], self.p['outcol'])


class Apply(BBS):
    """
    Action to apply solutions
    """
    def __init__(self, op_parset, vis_datamap, p, parmdb_datamap, prefix=None,
        direction=None, clean=True, index=None):
        super(Apply, self).__init__(op_parset, vis_datamap, p,
            parmdb_datamap=parmdb_datamap,
            prefix=prefix, direction=direction, clean=clean, index=index,
            name='Apply')


class Solve(BBS):
    """
    Action to solve for solutions
    """
    def __init__(self, op_parset, vis_datamap, p, model_datamap=None,
        parmdb_datamap=None, prefix=None, direction=None, clean=True,
        index=None):
        super(Solve, self).__init__(op_parset, vis_datamap, p,
            model_datamap=model_datamap, parmdb_datamap=parmdb_datamap,
            prefix=prefix, direction=direction, clean=clean, index=index,
            name='Solve')


class Subtract(BBS):
    """
    Action to subtract sources
    """
    def __init__(self, op_parset, vis_datamap, p, model_datamap=None,
        parmdb_datamap=None, prefix=None, direction=None, clean=True,
        index=None):
        super(Subtract, self).__init__(op_parset, vis_datamap, p,
            model_datamap=model_datamap, parmdb_datamap=parmdb_datamap,
            prefix=prefix, direction=direction, clean=clean, index=index,
            name='Subtract')
