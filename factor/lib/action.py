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
        factor_working_dir = op_parset['dir_working']

        # Set up directories needed by every action
        self.vis_dir = '{0}/visdata/'.format(factor_working_dir)

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

        self.pipeline_executable = '{0}/bin/genericpipeline.py'.format(
            self.op_parset['lofarroot'])

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
        Writes the pipeline configuration parset (specific to given action)
        """
        raise NotImplementedError


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


