"""
General operation library

Contains the master class for operations
"""
import os
import logging
import socket
from factor.lib.context import Timer
from factor.lib.scheduler_mp import Scheduler
from factor import _logging
from jinja2 import Environment, FileSystemLoader

DIR = os.path.dirname(os.path.abspath(__file__))
env_parset = Environment(loader=FileSystemLoader(os.path.join(DIR, '..', 'pipeline',
    'parsets')))
env_config = Environment(loader=FileSystemLoader(os.path.join(DIR, '..', 'pipeline')))


class Operation(object):
    """
    Generic operation class.
    """
    def __init__(self, parset, bands, direction, name=None):
        """
        Create Operation object

        Parameters
        ----------
        parset : str
            Parset of operation
        bands : list of Band objects
            Bands for this operation
        direction : Direction object
            Direction for this operation
        name : str, optional
            Name of the action
        """
        self.parset = parset.copy()
        self.bands = bands
        self.name = name.lower()
        self.parset['op_name'] = name
        self.direction = direction
        _logging.set_level(self.parset['logging_level'])
        self.log = logging.getLogger('factor.{0}'.format(self.name))
        self.hostname = socket.gethostname()

        self.factor_working_dir = parset['dir_working']
        self.statebasename = '{0}/state/{1}-{2}'.format(self.factor_working_dir,
            self.name, self.direction.name)

        self.mapfile_dir = '{0}/datamaps/{1}/{2}'.format(self.factor_working_dir,
            self.name, self.direction.name)
        if not os.path.exists(self.mapfile_dir):
            os.makedirs(self.mapfile_dir)

        self.pipeline_run_dir = '{0}/pipeline/{1}/{2}'.format(self.factor_working_dir,
            self.name, self.direction.name)
        if not os.path.exists(self.pipeline_run_dir):
            os.makedirs(self.pipeline_run_dir)

        self.log_dir = '{0}/logs/{1}/{2}/'.format(self.factor_working_dir,
            self.name, self.direction.name)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)
        self.logbasename = os.path.join(self.log_dir, '{0}_{1}'.format(
            self.name, self.direction.name))

        self.factor_root_dir = os.path.join(DIR, '..')
        self.factor_pipeline_dir = os.path.join(self.factor_root_dir, 'pipeline')
        self.factor_parset_dir = os.path.join(self.factor_root_dir, 'parsets')
        self.factor_skymodel_dir = os.path.join(self.factor_root_dir, 'skymodels')
        self.pipeline_parset_template = env_parset.get_template('{0}_pipeline.parset'.
            format(self.name))
        self.pipeline_parset_file = os.path.join(self.pipeline_run_dir,
            'pipeline.parset')
        self.pipeline_config_template = env_config.get_template('pipeline.cfg')
        self.pipeline_config_file = os.path.join(self.pipeline_run_dir,
            'pipeline.cfg')

        self.cfg_dict = {'lofarroot': parset['lofarroot'],
                         'pythonpath': parset['lofarpythonpath'],
                         'factorroot': self.factor_root_dir,
                         'working_dir': self.factor_working_dir,
                         'runtime_dir': self.pipeline_run_dir}


    def write_mapfile(self, data_list, prefix=None, direction=None, band=None,
        index=None, host_list=None):
        """
        Write an operation datamap.

        Parameters
        ----------
        data_list : list of str
            List of files for datamap
        prefix : str, optional
            A prefix for the name
        direction : Direction object, optional
            A direction
        band : Band object, optional
            A band
        index : int, optional
            An index for the datamap
        host_list : list of str, optional
            List of hosts for datamap

        """
        from factor.lib.datamap_lib import write_mapfile

        if host_list is None:
            host_list = self.parset['cluster_specific']['node_list']

        mapfile = write_mapfile(data_list, self.name, prefix=prefix,
                direction=direction, band=band, index=index, host_list=host_list,
                working_dir=self.parset['dir_working'])

        return mapfile


    def setup(self):
        """
        Set up this operation
        """

        tmp = self.pipeline_parset_template.render(self.parms_dict)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)
        tmp = self.pipeline_config_template.render(self.cfg_dict)
        with open(self.pipeline_config_file, 'w') as f:
            f.write(tmp)


    def finalize(self):
        """
        Finalize this operation
        """
        pass


    def check_completed(self):
        """
        Checks whether operation has been run successfully before

        Returns
        -------
        all_done : bool
            True if all objects were successfully run

        """
        self.direction.load_state()
        if self.name in self.direction.completed_operations:
            all_done = True
        else:
            all_done = False

        return all_done


    def set_completed(self):
        """
        Sets the state for the operation objects (bands or directions)
        """
        self.direction.completed_operations.append(self.name)
        self.direction.save_state()
