"""
General action library

Contains the master class for actions. Actions handle the setup and running of
the pipeline.
"""

import logging
import subprocess
import os

class Action(object):
    """
    Generic action class
    """
    def __init__(self, op_parset, name):
        """
        Create Action object

        Parameters
        ----------
        op_parset : str
            Parset of calling operation
        name : str
            Name of the action
        """
        self.op_name = op_parset['op_name']
        self.name = name
        self.op_parset = op_parset.copy()
        self.pipeline_parset_file = None
        self.pipeline_config_file = None

        self.image_dir = 'images/{0}/'.format(self.op_name)
        if not os.path.exists(self.image_dir):
            os.makedirs(self.image_dir)

        self.vis_dir = 'visdata/'

        self.datamap_dir = 'datamaps/{0}/{1}/'.format(self.op_name, self.name)
        if not os.path.exists(self.datamap_dir):
            os.makedirs(self.datamap_dir)

        self.model_dir = 'models/{0}/'.format(self.op_name)
        if not os.path.exists(self.model_dir):
            os.makedirs(self.model_dir)

        self.parset_dir = 'parsets/{0}/{1}/'.format(self.op_name, self.name)
        if not os.path.exists(self.parset_dir):
            os.makedirs(self.parset_dir)

        self.pipeline_run_dir = 'pipeline/{0}/{1}/'.format(self.op_name, self.name)
        if not os.path.exists(self.pipeline_run_dir):
            os.makedirs(self.pipeline_run_dir)

        self.log = logging.getLogger('%s::%s' % (self.op_name, self.name))
        self.log_dir = 'logs/{0}/{1}/'.format(self.op_name, self.name)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        self.pipeline_executable = '{0}/bin/genericpipeline.py'.format(
            self.op_parset['lofarroot'])


    def make_datamaps(self):
        """
        Makes the output datamap
        """
        raise NotImplementedError


    def make_pipeline_control_parset(self):
        """
        Makes the pipeline control parset
        """
        raise NotImplementedError


    def make_pipeline_config_parset(self):
        """
        Makes the pipeline configuration parset
        """
        raise NotImplementedError


    def run(self):
        """
        Runs the generic pipeline with parset
        """
        cmd = 'python {0} {1} -d -c {2}'.format(self.pipeline_executable,
            self.pipeline_parset_file, self.pipeline_config_file)
        with open("{0}.out.log".format(self.logbasename), "wb") as out, \
            open("{0}.err.log".format(self.logbasename), "wb") as err:
            p = subprocess.Popen(cmd, shell=True, stdout=out, stderr=err)
            p.wait()


    def get_results(self):
        """
        Return results data map
        """
        return self.outupt_datamap


