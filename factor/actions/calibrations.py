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

    Input data maps
    ---------------
    vis_datamap : Datamap
        Map of MS files
    model_datamap : Datamap, optional
        Map of sky models
    parmdb_datamap: Datamap, optional
        Map of parmdb files

    Output data maps
    ----------------
    parmdb_datamap: Datamap
        Map of parmdb files

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
        self.working_dir = self.vis_dir + '{0}/{1}/'.format(self.op_name, self.name)
        if self.direction is not None:
            self.working_dir += '{0}/'.format(self.direction.name)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)
        self.model_dir += '{0}/{1}/'.format(self.op_name, self.name)
        if self.direction is not None:
            self.model_dir += '{0}/'.format(self.direction.name)
        if not os.path.exists(self.model_dir):
            os.makedirs(self.model_dir)

        # Define parset name
        self.parset_file = self.parsetbasename + '.parset'
        self.templatename = '{0}_{1}.parset.tpl'.format(prefix, self.name.lower())

        # Set up all required files
        self.setup()


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import read_mapfile
        import copy

        self.p['vis_datamap'] = self.vis_datamap
        vis_files, vis_hosts = read_mapfile(self.vis_datamap)

        if self.model_datamap is None:
            emptyskymodel_list = []
            for v in vis_files:
                sm = os.path.join(self.model_dir, v+'.emptyskymodel')
                if not os.path.exists(sm):
                    os.system('touch {0}'.format(sm))
                emptyskymodel_list.append(sm)
            self.model_datamap = self.write_mapfile(emptyskymodel_list,
                prefix=self.prefix+'_empty_skymodels', index=self.index,
                host_list=vis_hosts, direction=self.direction)
        self.p['skymodel_datamap'] = self.model_datamap

        if self.parmdb_datamap is not None:
            self.p['parmdb_datamap'] = self.parmdb_datamap

        local_parmdb_files = [os.path.join(v, 'instrument') for v in vis_files]
        self.local_parmdb_datamap = self.write_mapfile(local_parmdb_files,
        	prefix=self.prefix+'_local_parmdbs', index=self.index,
        	host_list=vis_hosts, direction=self.direction)

        # Copy the input parmdb (if any) to local parmdb
        if self.parmdb_datamap is not None:
            parmdb_files, _ = read_mapfile(self.parmdb_datamap)
            if not os.path.exists(parmdb_files[0]):
                self.output_parmdb_datamap = copy.deepcopy(self.parmdb_datamap)
                self.parmdb_datamap = None
            else:
                for inp, outp in zip(parmdb_files, local_parmdb_files):
                    if os.path.exists(outp):
                        os.system('rm -rf {0}'.format(outp))
                    os.system('cp -r {0} {1}'.format(inp, outp))
                self.output_parmdb_datamap = None
        else:
            self.output_parmdb_datamap = None


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        if 'ncpu' not in self.p:
            self.p['ncpu'] = self.max_cpu
        self.p['outputdir'] = self.working_dir
        self.p['lofarroot'] = self.op_parset['lofarroot']
        self.p['parset'] = self.parset_file
        if 'flags' not in self.p:
            self.p['flags'] = ''
        if self.op_parset['use_ftw']:
            self.p['sources'] = '@MODEL_DATA'
        else:
            self.p['sources'] = ''

        template = env.get_template('bbs.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        template = env.get_template(self.templatename)
        tmp = template.render(self.p)
        with open(self.parset_file, 'w') as f:
            f.write(tmp)


    def get_results(self):
        """
        Return results

        Returns
        -------
        parmdb_datamap : DataMap
            Output solutions parmdb data map

        """
        from factor.lib.datamap_lib import read_mapfile

        default_parmdb_files, default_parmdb_hosts = read_mapfile(
            self.local_parmdb_datamap)

        if self.output_parmdb_datamap is not None:
            output_parmdb_files, output_parmdb_hosts = read_mapfile(
                self.output_parmdb_datamap)
            for inp, outp in zip(default_parmdb_files, output_parmdb_files):
                os.system('cp -r {0} {1}'.format(inp, outp))
            parmdb_datamap = self.output_parmdb_datamap
        elif self.parmdb_datamap is not None:
            parmdb_datamap = self.parmdb_datamap
        else:
            parmdb_datamap = self.local_parmdb_datamap

        return parmdb_datamap


class DPPP(Action):
    """
    Action to run DPPP Gaincal

    Input data maps
    ---------------
    vis_datamap : Datamap
        Map of MS files
    model_datamap : Datamap, optional
        Map of sky models
    parmdb_datamap: Datamap, optional
        Map of parmdb files

    Output data maps
    ----------------
    parmdb_datamap: Datamap, optional
        Map of parmdb files

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

        model_files, model_hosts = read_mapfile(self.model_datamap)
        model_flags = read_mapfile_flags(self.model_datamap)
        vis_files, vis_hosts = read_mapfile(self.vis_datamap)
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
        	parmdb_datamap=parmdb_datamap, prefix=prefix,
        	direction=direction, clean=clean, index=index, name='Apply')


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
    def __init__(self, op_parset, vis_datamap, p, model_datamap,
        parmdb_datamap, prefix=None, direction=None, clean=True,
        index=None):
        super(Subtract, self).__init__(op_parset, vis_datamap, p,
            model_datamap=model_datamap, parmdb_datamap=parmdb_datamap,
            prefix=prefix, direction=direction, clean=clean, index=index,
            name='Subtract')
