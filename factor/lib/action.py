"""
General action library

Contains the master class for actions. Actions handle the setup and running of
the pipeline.
"""

import logging
import subprocess
import os
from factor.lib.action_lib import make_basename


class Action(object):
    """
    Generic action class
    """
    def __init__(self, op_parset, name, prefix=None, direction=None,
        index=None):
        """
        Create Action object

        Parameters
        ----------
        op_parset : str
            Parset of calling operation
        name : str
            Name of the action
        prefix : str, optional
            Prefix to use for names
        direction : Direction object, optional
            Direction for this action
        index : int, optional
            Index of action

        """
        self.op_name = op_parset['op_name']
        self.name = name
        self.op_parset = op_parset.copy()
        self.prefix = prefix
        self.direction = direction
        self.index = index
        self.max_cpu = self.op_parset['ncpu']

        # Set up directories needed by every action
        factor_working_dir = op_parset['dir_working']
        self.vis_dir = '{0}/visdata/'.format(factor_working_dir)
        self.image_dir = '{0}/images/'.format(factor_working_dir)
        self.model_dir = '{0}/models/'.format(factor_working_dir)
        self.parmdb_dir = '{0}/parmdbs/'.format(factor_working_dir)

        self.datamap_dir = '{0}/datamaps/{1}/{2}/'.format(factor_working_dir,
            self.op_name, self.name)
        if self.direction is not None:
            self.datamap_dir += '{0}/'.format(self.direction.name)
        if not os.path.exists(self.datamap_dir):
            os.makedirs(self.datamap_dir)

        self.parset_dir = '{0}/parsets/{1}/{2}/'.format(factor_working_dir,
            self.op_name, self.name)
        if self.direction is not None:
            self.parset_dir += '{0}/'.format(self.direction.name)
        if not os.path.exists(self.parset_dir):
            os.makedirs(self.parset_dir)

        self.pipeline_run_dir = '{0}/pipeline/{1}/{2}/'.format(factor_working_dir,
            self.op_name, self.name)
        if self.direction is not None:
            self.pipeline_run_dir += '{0}/'.format(self.direction.name)
        if not os.path.exists(self.pipeline_run_dir):
            os.makedirs(self.pipeline_run_dir)

        self.log = logging.getLogger('%s::%s' % (self.op_name, self.name))
        self.log_dir = '{0}/logs/{1}/{2}/'.format(factor_working_dir,
            self.op_name, self.name)
        if self.direction is not None:
            self.log_dir += '{0}/'.format(self.direction.name)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        self.pipeline_executable = os.path.join(self.op_parset['piperoot'], 'bin',
            'genericpipeline.py')

        # Set up parset and script names needed by every action
        self.parsetbasename = self.parset_dir + make_basename(prefix,
            direction, index)
        self.pipeline_parset_file = self.parsetbasename + 'pipe.parset'
        self.pipeline_config_file = self.parsetbasename + 'pipe.cfg'
        self.pipeline_run_dir += make_basename(prefix, direction, index)
        if not os.path.exists(self.pipeline_run_dir):
            os.makedirs(self.pipeline_run_dir)
        self.logbasename = self.log_dir + make_basename(prefix, direction,
            index)


    def make_datamaps(self):
        """
        Makes the output datamaps (specific to given action)
        """
        raise NotImplementedError


    def write_mapfile(self, data_list, prefix=None, direction=None, index=None,
        host_list=None, flag_list=None, use_abs_path=True):
        """
        Write operation datamap

        Parameters
        ----------
        data_list : list of str
            List of files
        prefix : str
            A prefix for the name
        direction : Direction object or str, optional
            A direction name
        index : int, optional
            An index for the datamap
        flag_list : list of bools, optional
            List of datamap flags
        use_abs_path : bool, optional
            If True, make sure all files use absolute paths
        host_list : list, optional
            List of datamap host names

        """
        from factor.lib.datamap_lib import write_mapfile

        if host_list is None:
            host_list = self.op_parset['node_list']

        mapfile = write_mapfile(data_list, self.op_name,
        	action_name=self.name, prefix=prefix, direction=direction,
        	index=index, host_list=host_list, use_abs_path=use_abs_path,
        	working_dir=self.op_parset['dir_working'])

        return mapfile


    def make_pipeline_control_parset(self):
        """
        Makes the pipeline control parset (specific to given action)
        """
        raise NotImplementedError


    def set_pipeline_parameters(self):
        """
        Sets various pipeline parameters
        """
        self.op_parset['runtime_dir'] = self.pipeline_run_dir
        self.op_parset['working_dir'] = self.working_dir


    def make_pipeline_config_parset(self):
        """
        Writes the pipeline configuration parset
        """
        from jinja2 import Environment, FileSystemLoader

        DIR = os.path.dirname(os.path.abspath(__file__))
        env = Environment(loader=FileSystemLoader(os.path.join(DIR,
            '../actions/templates')))

        if os.path.basename(self.op_parset['clusterdesc']) == 'local.clusterdesc':
            self.op_parset['remote'] = '[remote]\n'\
                + 'method = local\n'\
                + 'max_per_node = {0}\n'.format(self.op_parset['ncpu'])
        else:
            self.op_parset['remote'] = ''

        template = env.get_template('pipeline.cfg.tpl')
        tmp = template.render(self.op_parset)
        with open(self.pipeline_config_file, 'w') as f:
            f.write(tmp)


    def run(self):
        """
        Runs the pipeline
        """
        self.make_datamaps()
        self.set_pipeline_parameters()
        self.make_pipeline_control_parset()
        self.make_pipeline_config_parset()

        cmd = 'python {0} {1} -d -c {2}'.format(self.pipeline_executable,
            self.pipeline_parset_file, self.pipeline_config_file)
        with open("{0}.out.log".format(self.logbasename), "wb") as out, \
            open("{0}.err.log".format(self.logbasename), "wb") as err:
            p = subprocess.Popen(cmd, shell=True, stdout=out, stderr=err)
            p.wait()

        return self.get_results()


    def get_results(self):
        """
        Return results data map
        """
        return self.output_datamap


