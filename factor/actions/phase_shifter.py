from factor.lib.action import action
import subprocess

commands = [ 'sleep 0.9', 'sleep 1.1']

class phase_shifter(action):
    """
    Implment the phase shifter action
    """

    def __init__(self, ms, direction):
        super(phase_shifter, self).__init__(name = 'Phase shifter')

    def run(self):
        for command in commands:
            #p = subprocess.Popen(command, stderr=outputfile, stdout=outputfile, shell=True)
            p = subprocess.Popen(command, shell=True)
            p.wait()
