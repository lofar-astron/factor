"""
Action: Example
An example action, to be used as a template
"""

from factor.lib.action import action

class example_action(action):
    """
    Implment the phase shifter action
    """

    def __init__(self, op_name, ms, direction):
        super(phase_shifter, self).__init__(op_name, name = 'example_action')

    def run(self):
        # these commands are run sequentially!
        commands = [ 'command 1', 'command 2', 'command 3']
        for command in commands:
            self.exec_cmd(command)
