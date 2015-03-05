"""
General operation library

Contains the master class for operations
"""
import os
import logging
from factor.lib.context import Timer
from factor.lib.scheduler import Scheduler

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
        self.log = logging.getLogger(self.name)
        self.s = Scheduler(parset['ncpu'], name=name)

        factor_working_dir = parset['dir_working']
        if self.direction is not None:
            if type(self.direction) is list:
                self.statebasename = []
                for d in self.direction:
                    self.statebasename.append('{0}/state/{1}-{2}'.format(factor_working_dir,
                        self.name, d.name)            )
            else:
                self.statebasename = '{0}/state/{1}-{2}'.format(factor_working_dir,
                    self.name, self.direction.name)
        else:
            self.statebasename = '{0}/state/{1}'.format(factor_working_dir,
                self.name)
        self.mapbasename = '{0}/datamaps/{1}/'.format(factor_working_dir, self.name)
        if not os.path.exists(self.mapbasename):
            os.makedirs(self.mapbasename)


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
            self.log.info('<-- Operation %s started (direction(s): %s)' % (self.name, dirstr))


    def write_datamap(self, data_list, prefix=None, direction=None, index=None,
        host_list=None):
        """
        Write operation datamap
        """
        from factor.lib.datamap_lib import write_mapfile

        if host_list is None:
            host_list = self.parset['node_list']

        mapfile = write_mapfile(data_list, self.name, prefix=prefix,
                direction=direction, index=index, host_list=host_list,
                working_dir=self.parset['dir_working']))

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
        with Timer(self.log):
            self.setup()
            self.run_steps()
            self.finalize()


    def finalize(self):
        """
        Finalize the operation
        """
        # Set the operation completion state. How to determine if all steps
        # completed successfully? Check pipeline state for each step?
        self.set_state(0)
        if self.direction is None:
            self.log.info('--> Operation %s finished' % self.name)
        else:
            if type(self.direction) is list:
                dirstr = ', '.join([d.name for d in self.direction])
            else:
                dirstr = self.direction.name
            self.log.info('--> Operation %s finished (direction(s): %s)' % (self.name, dirstr))


    def set_state(self, returncode):
        """
        Set success or failure state
        """
        if type(self.statebasename) is list:
            for s in self.statebasename:
                success_file = s + '.done'
                failure_file = s + '.failed'
            if returncode == 0:
                state_file = success_file
                if os.path.exists(failure_file):
                    os.remove(failure_file)
            else:
                state_file = failure_file
                if os.path.exists(success_file):
                    os.remove(success_file)
            with open(state_file, 'w') as f:
                f.write(' ')
            if returncode != 0 and self.exit_on_error:
                import sys
                sys.exit(1)
        else:
            success_file = self.statebasename + '.done'
            failure_file = self.statebasename + '.failed'
            if returncode == 0:
                state_file = success_file
                if os.path.exists(failure_file):
                    os.remove(failure_file)
            else:
                state_file = failure_file
                if os.path.exists(success_file):
                    os.remove(success_file)
            with open(state_file, 'w') as f:
                f.write(' ')
            if returncode != 0 and self.exit_on_error:
                import sys
                sys.exit(1)


    def unset_state(self):
        """
        Unset any previous state
        """
        success_file = self.statebasename + '.done'
        failure_file = self.statebasename + '.failed'
        if os.path.exists(failure_file):
            os.remove(failure_file)
        if os.path.exists(success_file):
            os.remove(success_file)

