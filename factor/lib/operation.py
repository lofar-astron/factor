"""
General operation library

Contains the master class for all operations
"""
import os
import logging
import socket
import numpy as np
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

    An operation is simply a generic pipeline that performs a part of the facet
    calibration. The corresponding operation object holds the pipeline settings,
    populates the pipeline config and parset templates, and updates the direction
    object with variables needed by later operations.

    Parameters
    ----------
    parset : dict
        Parset of operation
    bands : list of Band objects
        Bands for this operation
    direction : Direction object
        Direction for this operation
    name : str, optional
        Name of the operation

    """
    def __init__(self, parset, bands, direction, name=None):
        self.parset = parset.copy()
        self.bands = bands
        self.name = name.lower()
        self.parset['op_name'] = name
        self.direction = direction
        _logging.set_level(self.parset['logging_level'])
        self.log = logging.getLogger('factor:{0}'.format(self.name))
        self.hostname = socket.gethostname()
        self.node_list = parset['cluster_specific']['node_list']

        # Working directory
        self.factor_working_dir = parset['dir_working']

        # Pipeline runtime and working dirs (pipeline makes subdir here with
        # name of direction)
        self.pipeline_runtime_dir = os.path.join(self.factor_working_dir, 'results',
            self.name)
        self.pipeline_working_dir = self.pipeline_runtime_dir
        create_directory(self.pipeline_runtime_dir)

        # Directory that holds the mapfiles
        self.pipeline_mapfile_dir = os.path.join(self.pipeline_runtime_dir,
            self.direction.name, 'mapfiles')
        create_directory(self.pipeline_mapfile_dir)

        # Directory in the runtime dir that holds parset and config files (also
        # the results of the pipeline)
        self.pipeline_parset_dir = os.path.join(self.pipeline_runtime_dir,
            self.direction.name)
        create_directory(self.pipeline_parset_dir)

        # Directory that holds the mapfiles
        self.pipeline_mapfile_dir = os.path.join(self.pipeline_runtime_dir,
            self.direction.name, 'mapfiles')
        create_directory(self.pipeline_mapfile_dir)

        # Local scratch directory and corresponding node recipes
        if self.parset['cluster_specific']['dir_local'] is None:
            # Not specified, specify scratch directory in normal work directory
            self.local_scratch_dir = os.path.join(self.pipeline_working_dir,
                self.direction.name)
            self.dppp_nodescript = 'executable_args'
        elif self.parset['cluster_specific']['clusterdesc_file'].lower() == 'pbs':
            # PBS = "system in Hamburg" -> use special NDPPP nodescript
            self.local_scratch_dir = self.parset['cluster_specific']['dir_local']
            self.dppp_nodescript = 'dppp_scratch'
        else:
            # other: use given scratch directory an standard nodescrit
            self.local_scratch_dir = self.parset['cluster_specific']['dir_local']
            self.dppp_nodescript = 'executable_args'

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
        self.cfg_dict = {'lofarroot': parset['cluster_specific']['lofarroot'],
                         'pythonpath': parset['cluster_specific']['lofarpythonpath'],
                         'factorroot': self.factor_root_dir,
                         'pipeline_working_dir': self.pipeline_working_dir,
                         'pipeline_runtime_dir': self.pipeline_runtime_dir,
                         'casa_executable': parset['casa_executable'],
                         'wsclean_executable': parset['wsclean_executable'],
                         'image2fits_executable': parset['image2fits_executable'],
                         'dppp_nodescript': self.dppp_nodescript}

        # Define global parameters needed by all pipeline parsets. Other,
        # pipeline-specific, parameters should be defined in the subclasses by
        # updating this dictionary
        self.parms_dict = {'parset_dir': self.factor_parset_dir,
                           'skymodel_dir': self.factor_skymodel_dir,
                           'mapfile_dir': self.pipeline_mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'script_dir': self.factor_script_dir,
                           'local_dir': self.local_scratch_dir,
                           'pipeline_parset_dir': self.pipeline_parset_dir,
                           'hosts': self.node_list}

        # Add cluster-related info
        if self.parset['cluster_specific']['clustertype'] == 'local':
            self.cfg_dict['remote'] = '[remote]\n'\
                + 'method = local\n'\
                + 'max_per_node = {0}\n'.format(self.direction.max_cpus_per_node)
        elif self.parset['cluster_specific']['clustertype'] == 'juropa_slurm':
            self.cfg_dict['remote'] = '[remote]\n'\
                + 'method = slurm_srun\n'\
                + 'max_per_node = {0}\n'.format(self.direction.max_cpus_per_node)
        elif self.parset['cluster_specific']['clustertype'] == 'pbs':
            self.cfg_dict['remote'] = ''
        else:
            self.log.error('Could not determine the nature of your cluster!')
            sys.exit(1)

        # an absolute path in ...['clusterdesc'] will overrule the the "working_dir"
        self.cfg_dict['clusterdesc'] = os.path.join(self.factor_working_dir,
            self.parset['cluster_specific']['clusterdesc'])


    def update_dicts(self):
        self.cfg_dict.update(self.direction.__dict__)
        self.parms_dict.update(self.direction.__dict__)


    def setup(self):
        """
        Set up this operation

        This involves just filling the pipeline config and parset templates.
        Generally, this does not need to be re-defined in the subclasses
        unless the operation has non-standard template names
        """
        # Update the dictionaries with the attributes of the operation's
        # direction object. Any attributes set in the direction object that are
        # also in the parms_dict will be set to those of the direction object
        # (e.g., 'max_cpus_per_node', which is set in the direction object by
        # factor.cluster.divide_nodes() will override the value set above)
        self.update_dicts()

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


    def check_started(self):
        """
        Checks whether operation has been started (but not necessarily
        completed) before for this direction

        Returns
        -------
        all_done : bool
            True if operation was started on this direction

        """
        has_state = self.direction.load_state()
        if has_state:
            if self.name in self.direction.started_operations:
                return True
            else:
                return False
        else:
            return False


    def check_completed(self):
        """
        Checks whether operation has been run successfully before for this
        direction

        Returns
        -------
        all_done : bool
            True if operation was successfully run on this direction

        """
        has_state = self.direction.load_state()
        if has_state:
            if self.name in self.direction.completed_operations:
                return True
            else:
                return False
        else:
            return False


    def set_started(self):
        """
        Sets the started state for the operation
        """
        if self.name not in self.direction.started_operations:
            self.direction.started_operations.append(self.name)
            self.direction.save_state()


    def set_completed(self):
        """
        Sets the completed state for the operation
        """
        if self.name not in self.direction.completed_operations:
            self.direction.completed_operations.append(self.name)
            self.direction.save_state()
