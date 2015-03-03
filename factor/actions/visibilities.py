"""
Module that holds all visibility-related actions

Classes
-------
DPPP : Action
    Runs DPPP
Average : Action
    Averages visibilities
PhaseShift : Action
    Phase shifts visibilities
Concatenate : Action
    Concatenates visibilities

"""

from factor.lib.action import Action
from factor.lib.action_lib import make_basename
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class DPPP(Action):
    """
    Action to run DPPP

    Input data maps
    ---------------
    input_datamap : Datamap
        Map of vis data files

    Output data maps
    ----------------
    output_datamap : Datamap
        Map of output vis files

    """
    def __init__(self, op_parset, input_datamap, p, prefix=None,
        direction=None, clean=True, index=None, name='DPPP'):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map for action
        p : dict
            Input parset dict defining model and pipeline parameters
        prefix : str, optional
            Prefix to use for model names
        direction : Direction object, optional
            Direction for this model
        clean : bool, optional
            Remove unneeded files?
        index : int, optional
            Index of action

        """
        super(DPPP, self).__init__(op_parset, name=name, prefix=prefix,
            direction=direction, index=index)

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'run_dppp'
        self.clean = clean
        self.working_dir = self.vis_dir + '{0}/{1}/'.format(self.op_name, self.name)
        if self.direction is not None:
            self.working_dir += '{0}/'.format(self.direction.name)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        if self.name != 'Concatenate':
            self.p['input_datamap'] = self.input_datamap

            msnames = read_mapfile(self.input_datamap)
            if self.index is not None:
                indstr = '{0}'.format(self.index)
            else:
                indstr = ''
            output_files = [self.working_dir + os.path.splitext(os.path.basename(ms))[0]
                + '_{0}{1}.ms'.format(self.name.lower(), indstr) for ms in msnames]
        else:
            # Handle concatenation separately, as we need a data map with a
            # single "file" (and hence must disable use_abs_path below).
            msnames = read_mapfile(self.input_datamap)
            concat_file_list = ['[' + ','.join([msname for msname in msnames]) + ']']
            self.p['input_datamap'] = write_mapfile(concat_file_list,
                self.op_name, self.name, prefix=self.prefix+'_input',
                direction=self.direction, working_dir=self.op_parset['dir_working'],
                use_abs_path=False)

            ms = msnames[0]
            output_files = [self.working_dir + os.path.splitext(os.path.basename(ms))[0]
                + '_{0}.ms'.format(self.name.lower())]

        self.p['output_datamap'] = write_mapfile(output_files,
            self.op_name, self.name, prefix=self.prefix+'_output',
            direction=self.direction, working_dir=self.op_parset['dir_working'])


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        self.p['lofarroot'] = self.op_parset['lofarroot']

        template = env.get_template('{0}.pipeline.parset.tpl'.format(self.name.lower()))
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)


    def make_pipeline_config_parset(self):
        """
        Writes the pipeline configuration parset
        """
        template = env.get_template('pipeline.cfg.tpl')
        tmp = template.render(self.op_parset)
        with open(self.pipeline_config_file, 'w') as f:
            f.write(tmp)


    def get_results(self):
        """
        Return results
        """
        return self.p['output_datamap']


class Average(DPPP):
    """
    Action to average visibilities
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        super(Average, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Average')


class PhaseShift(DPPP):
    """
    Action to phase shift visibilities
    """
    def __init__(self, op_parset, input_datamap, p, direction, prefix=None,
        clean=True, index=None):
        super(PhaseShift, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='PhaseShift')
        self.p['ra'] = direction.ra
        self.p['dec'] = direction.dec


class Concatenate(DPPP):
    """
    Action to concatenate visibilities
    """
    def __init__(self, op_parset, input_datamap, p, direction, prefix=None,
        clean=True, index=None):
        super(Concatenate, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Concatenate')

