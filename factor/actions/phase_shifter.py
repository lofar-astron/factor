from factor.lib.action import action
import subprocess

commands = [ 'echo aaa', 'sleep 0.1', 'glfgsdfsd']

class phase_shifter(action):
    """
    Implment the phase shifter action
    """

    def __init__(self, ms, direction):
        super(phase_shifter, self).__init__(name = 'phase_shifter')

    def run(self):
        for command in commands:
            self.exec_cmd(command)
