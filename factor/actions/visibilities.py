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
Split : Action
    Splits visibilities

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
        direction=None, band=None, clean=True, index=None, name='DPPP'):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map for action
        p : dict
            Input parset dict defining dppp and pipeline parameters
        prefix : str, optional
            Prefix to use
        direction : Direction object, optional
            Direction for this action
        band : Band object, optional
            Band for this action
        clean : bool, optional
            Remove unneeded files?
        index : int, optional
            Index of action

        """
        super(DPPP, self).__init__(op_parset, name=name, prefix=prefix,
            direction=direction, band=band, index=index)

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'run_dppp'
        self.clean = clean
        self.working_dir = self.vis_dir + '{0}/{1}/'.format(self.op_name, self.name)
        if self.direction is not None:
            self.working_dir += '{0}/'.format(self.direction.name)
        if self.band is not None:
            self.working_dir += '{0}/'.format(self.band.name)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        # Set up all required files
        self.setup()


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import read_mapfile

        msnames, hosts = read_mapfile(self.input_datamap)
        if self.name != 'Concatenate':
            self.p['input_datamap'] = self.input_datamap
            if self.index is not None:
                indstr = '{0}'.format(self.index)
            else:
                indstr = ''
            output_files = [self.working_dir + os.path.splitext(os.path.basename(ms))[0]
                + '_{0}{1}.ms'.format(self.name.lower(), indstr) for ms in msnames]
        else:
            # Handle concatenation separately, as we need a data map with a
            # single "file" (and hence must disable use_abs_path below).
            concat_file_list = ['[' + ','.join([msname for msname in msnames]) + ']']
            hosts = [hosts[0]]
            self.p['input_datamap'] = self.write_mapfile(concat_file_list,
                prefix=self.prefix+'_input', direction=self.direction, band=self.band,
                use_abs_path=False, index=self.index, host_list=hosts)
            ms = msnames[0]
            output_files = [self.working_dir + os.path.splitext(os.path.basename(ms))[0]
                + '_{0}.ms'.format(self.name.lower())]

        self.p['output_datamap'] = self.write_mapfile(output_files,
            prefix=self.prefix+'_output', direction=self.direction,
            index=self.index, band=self.band, host_list=hosts)


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        if 'ncpu' not in self.p:
            self.p['ncpu'] = self.max_cpu
        if 'n_per_node' not in self.p:
            self.p['n_per_node'] = self.max_cpu
        self.p['lofarroot'] = self.op_parset['lofarroot']

        template = env.get_template('{0}.pipeline.parset.tpl'.format(self.name.lower()))
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
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
        band=None, clean=True, index=None):
        super(Average, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, band=band, clean=clean, index=index, name='Average')


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

        # Set up all required files
        self.setup()


class Concatenate(DPPP):
    """
    Action to concatenate visibilities
    """
    def __init__(self, op_parset, input_datamap, p, direction=None, prefix=None,
        clean=True, index=None):
        super(Concatenate, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, clean=clean, index=index, name='Concatenate')


class Split(DPPP):
    """
    Action to split off visibilities
    """
    def __init__(self, op_parset, input_datamap, p, direction=None, prefix=None,
        band=None, clean=True, index=None):
        super(Split, self).__init__(op_parset, input_datamap, p, prefix=prefix,
            direction=direction, band=band, clean=clean, index=index, name='Split')

