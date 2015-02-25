"""
Module that holds all visibility-related actions

Classes
-------
DPPP : Action
    Runs DPPP
Average : Action
    Averages visibilities
PhaseShift : Action
    Phase shifts visibilities
Concatenate : Action
    Concatenates visibilities

"""

from factor.lib.action import Action
from factor.lib.action_lib import make_basename
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class DPPP(Action):
    """
    Action to run DPPP
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None,
        direction=None, clean=True, index=None, name='DPPP'):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map for action
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
        super(DPPP, self).__init__(op_parset, name=name, prefix=prefix,
            direction=direction, index=index)

        # Store input parameters
        self.input_datamap = vis_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'run_dppp'
        self.clean = clean
        factor_working_dir = op_parset['dir_working']
        if self.direction is None:
            self.working_dir = '{0}/visdata/{1}/'.format(factor_working_dir,
                self.op_name)
        else:
            self.working_dir = '{0}/visdata/{1}/{2}/'.format(factor_working_dir,
                self.op_name, self.direction)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        self.p['input_datamap'] = self.input_datamap

        msnames = read_mapfile(self.input_datamap)
        output_files = [self.working_dir + os.path.splitext(os.path.basename(ms))[0]
            + '_avg.ms' for ms in msnames]

        self.p['output_datamap'] = write_mapfile(output_files,
            self.op_name, self.name, prefix=self.prefix+'_output',
            direction=self.direction, working_dir=self.op_parset['dir_working'])


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        self.p['lofarroot'] = self.op_parset['lofarroot']

        template = env.get_template('{0}.pipeline.parset.tpl'.format(self.name.lower()))
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
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
        return


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



class Average(DPPP):
    """
    Action to average visibilities
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        super(Average, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Average')


class PhaseShift(DPPP):
    """
    Action to phase shift visibilities
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        super(PhaseShift, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='PhaseShift')


class Concatenate(DPPP):
    """
    Action to concatenate visibilities
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        super(Concat, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Concat')

