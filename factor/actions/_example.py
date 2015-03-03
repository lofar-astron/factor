"""
Module that holds an example actions

Action: Example
An example action, to be used as a template
"""

from factor.lib.action import Action

class ExampleAction(Action):
    """
    An example action
    """

    def __init__(self, op_parset, input_datamap):
        """
        Create Action object

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map of file(s) to process
        """
        super(ExampleAction, self).__init__(op_parset, name='ExampleAction')
        self.input_datamap = input_datamap

        # Make data maps
        self.make_datamaps()

        # Make pipeline parsets
        self.set_pipeline_parameters()
        self.make_pipeline_control_parset()
        self.make_pipeline_parset()

        # Run the pipeline
        self.run()
