"""
General operation library

Contains the master class for operations
"""
import os
import logging
from factor.lib.context import Timer


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
        factor_working_dir = parset['dir_working']

        if self.direction is not None:
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
            self.log.info('<-- Operation %s started (direction: %s)' % (self.name, self.direction.name))


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
            self.log.info('--> Operation %s finished (direction: %s)' % (self.name, self.direction.name))


    def set_state(self, returncode):
        """
        Set success or failure state
        """
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

