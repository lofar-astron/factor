from factor.lib.action import action

class phase_shifter(action):
    """
    Implment the phase shifter action
    """

    def __init__(self, op_name, ms, direction):
        super(phase_shifter, self).__init__(op_name, name = 'phase_shifter')

    def run(self):
        commands = [ 'echo aaa', 'sleep 1', 'glfgsdfsd']
        self.log.debug('Running phase shifter')
        for command in commands:
            self.exec_cmd(command)
