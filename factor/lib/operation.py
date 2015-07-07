"""
General operation library

Contains the master class for all operations
"""
import os
import logging
import socket
from factor.lib.context import Timer
from factor.lib.scheduler_mp import Scheduler
from factor import _logging
from jinja2 import Environment, FileSystemLoader
from lofarpipe.support.utilities import create_directory

DIR = os.path.dirname(os.path.abspath(__file__))
env_parset = Environment(loader=FileSystemLoader(os.path.join(DIR, '..', 'pipeline',
    'parsets')))
env_config = Environment(loader=FileSystemLoader(os.path.join(DIR, '..', 'pipeline')))


class Operation(object):
    """
    Generic operation class
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
            Name of the operation
        """
        self.parset = parset.copy()
        self.bands = bands
        self.name = name.lower()
        self.parset['op_name'] = name
        self.direction = direction
        _logging.set_level(self.parset['logging_level'])
        self.log = logging.getLogger('factor.{0}'.format(self.name))
        self.hostname = socket.gethostname()
        self.node_list = parset['cluster_specific']['node_list']
        self.max_cpus_per_node = parset['cluster_specific']['ncpu']
        self.max_percent_memory = parset['cluster_specific']['fmem'] * 100.0
        self.ndir_per_node = parset['cluster_specific']['ndir_per_node']

        # Below are paths for output directories
        self.factor_working_dir = parset['dir_working']
        # Base name of state file
        self.statebasename = os.path.join(self.factor_working_dir,
            'state', '{0}-{1}'.format(self.name, self.direction.name))
        # Directory that holds important mapfiles in a convenient place
        self.mapfile_dir = os.path.join(self.factor_working_dir, 'datamaps',
            self.name, self.direction.name)
        create_directory(self.mapfile_dir)
        # Pipeline runtime dir (pipeline makes subdir here with name of direction)
        self.pipeline_runtime_dir = os.path.join(self.factor_working_dir, 'results',
            self.name)
        create_directory(self.pipeline_runtime_dir)
        # Directory that holds parset and config files
        self.pipeline_parset_dir = os.path.join(self.pipeline_runtime_dir,
            self.direction.name)
        create_directory(self.pipeline_parset_dir)
        # Pipeline working dir (pipeline makes subdir here with name of direction)
        self.pipeline_working_dir = os.path.join(self.factor_working_dir, 'results',
            self.name)
        create_directory(self.pipeline_working_dir)
        # Directory that holds logs in a convenient place
        self.log_dir = os.path.join(self.factor_working_dir, 'logs', self.name)
        create_directory(self.log_dir)
        # Log name used for logs in log_dir
        self.logbasename = os.path.join(self.log_dir, self.direction.name)

        # Below are paths for scripts, etc. in the Factor install directory
        self.factor_root_dir = os.path.split(DIR)[0]
        self.factor_pipeline_dir = os.path.join(self.factor_root_dir, 'pipeline')
        self.factor_script_dir = os.path.join(self.factor_root_dir, 'scripts')
        self.factor_parset_dir = os.path.join(self.factor_root_dir, 'parsets')
        self.factor_skymodel_dir = os.path.join(self.factor_root_dir, 'skymodels')

        # Below are the templates and output paths for the pipeline parset and
        # config files. These may need to be re-defined in the subclasses
        # if the operation has non-standard template names
        self.pipeline_parset_template = '{0}_pipeline.parset'.format(self.name)
        self.pipeline_parset_file = os.path.join(self.pipeline_parset_dir,
            'pipeline.parset')
        self.pipeline_config_template = 'pipeline.cfg'
        self.pipeline_config_file = os.path.join(self.pipeline_parset_dir,
            'pipeline.cfg')

        # Define parameters needed for the pipeline config.
        self.cfg_dict = {'lofarroot': parset['lofarroot'],
                         'pythonpath': parset['lofarpythonpath'],
                         'factorroot': self.factor_root_dir,
                         'genericpiperoot': os.path.dirname(parset['genericpipeline_executable']).split('/bin')[0],
                         'pipeline_working_dir': self.pipeline_working_dir,
                         'pipeline_runtime_dir': self.pipeline_runtime_dir,
                         'max_cpus_per_node': self.max_cpus_per_node,
                         'casa_executable': parset['casa_executable'],
                         'wsclean_executable': parset['wsclean_executable'],
#                          'chgcentre_executable': parset['chgcentre_executable'],
                         'losoto_executable': parset['losoto_executable'],
                         'H5parm_importer_executable': parset['H5parm_importer_executable'],
                         'H5parm_exporter_executable': parset['H5parm_exporter_executable']}

        # Define global parameters needed by all pipeline parsets. Other,
        # pipeline-specific, parameters should be defined in the subclasses by
        # updating this dictionary
        self.parms_dict = {'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'script_dir': self.factor_script_dir,
                           'hosts': self.node_list,
                           'max_cpus_per_node': self.max_cpus_per_node,
                           'max_percent_memory' : self.max_percent_memory}
        self.parms_dict.update(self.direction.__dict__) # add all direction attributes

        # Add cluster-related info
        if os.path.basename(self.parset['cluster_specific']['clusterdesc']) == 'local.clusterdesc':
            self.cfg_dict['remote'] = '[remote]\n'\
                + 'method = local\n'\
                + 'max_per_node = {0}\n'.format(self.max_cpus_per_node)
        else:
            self.cfg_dict['remote'] = ''
        self.cfg_dict['clusterdesc'] = os.path.join(self.factor_working_dir,
            self.parset['cluster_specific']['clusterdesc'])


    def setup(self):
        """
        Set up this operation

        This involves just filling the pipeline config and parset templates.
        Generally, this does not need to be re-defined in the subclasses
        unless the operation has non-standard template names
        """
        self.pipeline_parset_template = env_parset.get_template(self.pipeline_parset_template)
        tmp = self.pipeline_parset_template.render(self.parms_dict)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        self.pipeline_config_template = env_config.get_template(self.pipeline_config_template)
        tmp = self.pipeline_config_template.render(self.cfg_dict)
        with open(self.pipeline_config_file, 'w') as f:
            f.write(tmp)


    def finalize(self):
        """
        Finalize this operation

        This should be defined in the subclasses if needed
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
            return True
        else:
            return False


    def set_completed(self):
        """
        Sets the state for the operation
        """
        self.direction.completed_operations.append(self.name)
        self.direction.save_state()
