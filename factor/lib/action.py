"""
General action library, contains the master class for actions
Commands contained in an action are run in sequential mode.
"""

import logging
import subprocess

class action( object ):
    """
    Generic action class
    """
    def __init__(self, op_name = None, name = None):
        self.op_name = op_name
        self.name = name
        self.log = logging.getLogger('%s::%s' % (self.op_name, self.name))

    def run(self):
        raise NotImplementedError

    def get_results(self):
        raise NotImplementedError

    def exec_cmd(self, cmd):
        # TODO: find a better name for logs, so they don't clash
        with open("log/%s-%s.out.log" % (self.name, cmd.partition(' ')[0]),"wb") as out, \
             open("log/%s-%s.err.log" % (self.name, cmd.partition(' ')[0]),"wb") as err:
            p = subprocess.Popen(cmd, shell=True, stdout=out, stderr=err)
            p.wait()

