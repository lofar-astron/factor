"""
Module that holds all solution-related actions

Classes
-------
Losoto : Action
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


class Losoto(Action):
    """
    Action to run LoSoTo
    """
    def __init__(self, op_parset, vis_datamap, p, parmdb_datamap, prefix=None,
        direction=None, clean=True, index=None, name='Losoto'):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        vis_datamap : data map
            Input data map for vis files
        p : dict
            Input parset dict defining model and pipeline parameters
        parmdb_datamap : data map
            Input data map for parmdb files
        prefix : str, optional
            Prefix to use for model names
        direction : Direction object, optional
            Direction for this model
        clean : bool, optional
            Remove unneeded files?
        index : int, optional
            Index of action

        """
        super(Losoto, self).__init__(op_parset, name=name, prefix=prefix,
            direction=direction, index=index)

        # Store input parameters
        self.vis_datamap = vis_datamap
        self.parmdb_datamap = parmdb_datamap
        self.p = p.copy()
        if prefix is None:
            prefix = 'run_losoto'
        self.clean = clean
        self.working_dir = self.parmdb_dir + '{0}/{1}/'.format(self.op_name, self.name)
        if self.direction is not None:
            self.working_dir += '{0}/'.format(self.direction.name)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        # Define parset name
        self.parset_file = self.parsetbasename + '.parset'
        self.templatename = '{0}_{1}.parset.tpl'.format(prefix, self.name.lower())

        # Copy parmdbs to instrument directory inside MS files
        parmdb_names = read_mapfile(self.parmdb_datamap)
        ms_names = read_mapfile(self.vis_datamap)
        for pn, mn in zip(parmdb_names, ms_names):
            if os.path.exists('{0}/instrument'.format(mn)):
                os.system('rm -rf {0}/instrument'.format(mn))
            os.system('cp -r {0} {1}/instrument'.format(pn, mn))


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        self.p['input_datamap'] = self.vis_datamap

        parmdb_names = read_mapfile(self.parmdb_datamap)
        h5parm_files = [parmdb_name + '.h5parm' for parmdb_name in parmdb_names]
        self.p['h5parm_datamap'] = write_mapfile(h5parm_files,
            self.op_name, self.name, prefix=self.prefix+'_h5parm',
            direction=self.direction, working_dir=self.op_parset['dir_working'])


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        self.p['lofarroot'] = self.op_parset['lofarroot']
        self.p['parset'] = self.parset_file

        template = env.get_template('{0}.pipeline.parset.tpl'.format(self.name.lower()))
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
        # Copy parmdbs from instrument directory inside MS files to final files
        parmdb_names = read_mapfile(self.parmdb_datamap)
        ms_names = read_mapfile(self.vis_datamap)
        for pn, mn in zip(parmdb_names, ms_names):
            if os.path.exists(pn):
                os.system('rm -rf {0}'.format(pn))
            os.system('cp -r {0}/{1}_instrument {2}'.format(mn, self.p['solset'], pn))

        return self.parmdb_datamap



class Smooth(Losoto):
    """
    Action to smooth and normalize solutions
    """
    def __init__(self, op_parset, vis_datamap, p, parmdb_datamap, prefix=None, direction=None,
        clean=True, index=None):
        super(Smooth, self).__init__(op_parset, vis_datamap, p, parmdb_datamap,
            prefix=prefix, direction=direction, clean=clean, index=index,
            name='Smooth')


class ResetPhases(Losoto):
    """
    Action to reset phases to zero
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        super(ResetPhases, self).__init__(op_parset, vis_datamap, p, parmdb_datamap,
            prefix=prefix, direction=direction, clean=clean, index=index,
            name='ResetPhases')
