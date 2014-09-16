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
    def __init__(self, name = None):
        self.name = name

    def run(self):
        raise NotImplementedError

    def exec_cmd(self, cmd):
        with open("log/%s-%s.out.log" % (self.name, cmd.partition(' ')[0]),"wb") as out, \
             open("log/%s-%s.err.log" % (self.name, cmd.partition(' ')[0]),"wb") as err:
            p = subprocess.Popen(cmd, shell=True, stdout=out, stderr=err)
            p.wait()

