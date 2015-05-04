"""
General operation library

Contains the master class for operations
"""
import os
import logging
import socket
from factor.lib.context import Timer
from factor.lib.scheduler import Scheduler
from factor import _logging

class Operation(object):
    """
    Generic operation class.

    All operations should be in a separate module. Every module must have a
    class called in the same way of the module which inherits from this class.
    """
    def __init__(self, parset, bands, direction=None, reset=False, name=None):
        """
        Create Operation object

        Parameters
        ----------
        parset : str
            Parset of operation
        bands : list of Band objects
            Bands for this operation
        direction : Direction object, optional
            Direction for this operation
        reset : bool, optional
            If True, reset the state of this operation
        name : str, optional
            Name of the action
        """
        self.parset = parset.copy()
        self.bands = bands
        self.name = name
        self.parset['op_name'] = name
        self.direction = direction
        self.reset = reset
        self.exit_on_error = True
        _logging.set_level(self.parset['logging_level'])
        self.log = logging.getLogger(self.name)
        self.s = Scheduler(parset['cluster_specific']['ncpu'], name=name,
            op_parset=self.parset)
        self.hostname = socket.gethostname()

        factor_working_dir = parset['dir_working']
        if self.direction is not None:
            if type(self.direction) is list:
                self.statebasename = []
                for d in self.direction:
                    self.statebasename.append('{0}/state/{1}-{2}'.format(
                        factor_working_dir, self.name, d.name))
            else:
                self.statebasename = '{0}/state/{1}-{2}'.format(factor_working_dir,
                    self.name, self.direction.name)
        else:
            self.statebasename = '{0}/state/{1}'.format(factor_working_dir,
                self.name)
        self.mapbasename = '{0}/datamaps/{1}/'.format(factor_working_dir, self.name)
        if not os.path.exists(self.mapbasename):
            os.makedirs(self.mapbasename)
        self.visbasename = '{0}/visdata/{1}/'.format(factor_working_dir, self.name)
        if not os.path.exists(self.visbasename):
            os.makedirs(self.visbasename)


    def setup(self):
        """
        Set up the operation
        """
        if self.direction is None:
            self.log.info('<-- Operation %s started' % self.name)
        else:
            if type(self.direction) is list:
                dirstr = ', '.join([d.name for d in self.direction])
            else:
                dirstr = self.direction.name
            self.log.info('<-- Operation %s started (direction(s): %s)' %
                (self.name, dirstr))


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


    def run_steps(self):
        """
        Define the operation's steps
        """
        raise(NotImplementedError)


    def run(self):
        """
        Run the operation
        """
        with Timer(self.log, 'operation'):
            self.setup()
            self.run_steps()
            self.finalize()


    def finalize(self):
        """
        Finalize the operation
        """
        # Set the operation completion state
        if self.direction is None:
            self.log.info('--> Operation %s finished' % self.name)
        else:
            if type(self.direction) is list:
                dirstr = ', '.join([d.name for d in self.direction])
            else:
                dirstr = self.direction.name
            self.log.info('--> Operation %s finished (direction(s): %s)' %
                (self.name, dirstr))


    def check_completed(self, obj_list):
        """
        Checks whether operation has been run successfully before

        Parameters
        ----------
        obj_list : list
            Band or Direction objects to check

        Returns
        -------
        all_done : bool
            True if all objects were successfully run

        """
        if type(obj_list) is not list:
            obj_list = [obj_list]

        for obj in obj_list:
            obj.load_state()
            if self.name in obj.completed_operations:
                all_done = True
            else:
                all_done = False
                break

        return all_done


    def set_completed(self, obj_list):
        """
        Sets the state for the operation objects

        Parameters
        ----------
        obj_list : list
            Band or Direction objects to check

        """
        if type(obj_list) is not list:
            obj_list = [obj_list]

        for obj in obj_list:
            obj.completed_operations.append(self.name)
            obj.save_state()
