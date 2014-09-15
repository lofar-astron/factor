"""
Action: Example
An example action, to be used as a template
"""

from factor.lib.action import action
import subprocess

# these commands are run sequentially!
commands = [ 'command 1', 'command 2', 'command 3']

class example_action(action):
    """
    Implment the phase shifter action
    """

    def __init__(self, ms, direction):
        super(phase_shifter, self).__init__(name = 'example_action')

    def run(self):
        for command in commands:
            p = subprocess.Popen(command, shell=True)
            p.wait()
