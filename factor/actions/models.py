"""
Module that holds all model-related actions

Classes
-------
MakeSkymodelFromModelImage : Action
    Makes a sky model from a CASA model image
MakeFacetSkymodel : Action
    Makes a sky model for a facet
MergeSkymodels : Action
    Merges two sky models
FFT : Action
    FFTs a model into an MS

"""

from factor.lib.action import Action
from factor.lib.action_lib import make_image_basename
from jinja2 import Environment, FileSystemLoader
import os

DIR = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(DIR, 'templates')))


class MakeSkymodelFromModelImage(Action):
    """
    Action to make a sky model from a CASA model image
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map for CASA model image
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
        super(MakeSkymodelFromModelImage, self).__init__(op_parset,
            'MakeSkymodelFromModelImage', prefix=prefix, direction=direction,
            index=index)

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'make_skymodel'
        self.clean = clean
        factor_working_dir = op_parset['dir_working']
        if self.direction is None:
            self.working_dir = '{0}/models/{1}/'.format(factor_working_dir,
                self.op_name)
        else:
            self.working_dir = '{0}/models/{1}/{2}/'.format(factor_working_dir,
                self.op_name, self.direction)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        # Define script name
        self.script_file = self.parsetbasename + 'make_skymodel_from_model_image.py'

        # Set up names for output data map
        modelbasenames = make_image_basename(self.input_datamap,
            direction=self.direction, prefix=self.prefix)
        self.modelbasenames = [self.working_dir+bn for bn in modelbasenames]


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        # Input is list of image basenames
        # Output is sky model files
        imagebasenames = read_mapfile(self.input_datamap)
        input_files = [bn+'.model' for bn in imagebasenames]
        output_files = [bn+'.skymodel' for bn in self.modelbasenames]

        self.p['input_datamap'] = write_mapfile(input_files,
            self.op_name, self.name, prefix=self.prefix+'_input',
            direction=self.direction, working_dir=self.op_parset['dir_working'])

        self.p['output_datamap'] = write_mapfile(output_files,
            self.op_name, self.name, prefix=self.prefix+'_output',
            direction=self.direction, working_dir=self.op_parset['dir_working'])


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        self.p['scriptname'] = os.path.abspath(self.script_file)
        self.p['outputdir'] = os.path.abspath(self.working_dir)
        template = env.get_template('make_skymodel_from_model_image.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        template = env.get_template('make_skymodel_from_model_image.tpl')
        tmp = template.render(self.p)
        with open(self.script_file, 'w') as f:
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
        Return skymodel names
        """
        return self.p['output_datamap']


class MakeFacetSkymodel(Action):
    """
    Action to make a sky model for a facet
    """
    def __init__(self, op_parset, input_datamap, p, direction, prefix=None,
        clean=True, index=None, cal_only=True):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map for full sky models
        p : dict
            Input parset dict defining model and pipeline parameters
        direction : Direction object
            Direction for this model
        prefix : str, optional
            Prefix to use for model names
        clean : bool, optional
            Remove unneeded files?
        index : int, optional
            Index of action
        cal_only : bool, optional
            If True, return sky model for calibrator only. If False, return
            sky model for entire facet

        """
        super(MakeFacetSkymodel, self).__init__(op_parset, 'MakeFacetSkymodel',
            prefix=prefix, direction=direction, index=index)

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'make_facet_skymodel'
        self.clean = clean
        self.cal_only = cal_only
        factor_working_dir = op_parset['dir_working']
        self.working_dir = '{0}/models/{1}/{2}/'.format(factor_working_dir,
            self.op_name, self.direction.name)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        # Define script name
        self.script_file = self.parsetbasename + 'make_facet_skymodels.py'

        # Set up names for output data map
        modelbasenames = make_image_basename(self.input_datamap,
            direction=self.direction, prefix=self.prefix)
        self.modelbasenames = [self.working_dir+bn for bn in modelbasenames]


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        self.p['input_datamap'] = self.input_datamap

        output_files = [bn + '_facet.skymodel' for bn in self.modelbasenames]
        self.p['output_datamap'] = write_mapfile(output_files,
            self.op_name, self.name, prefix=self.prefix+'_output',
            direction=self.direction, working_dir=self.op_parset['dir_working'])


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        self.p['scriptname'] = os.path.abspath(self.script_file)
        self.p['cal_only'] = self.cal_only
        self.p['vertices'] = self.direction.vertices
        template = env.get_template('make_facet_skymodel.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        template = env.get_template('make_facet_skymodel.tpl')
        tmp = template.render(self.p)
        with open(self.script_file, 'w') as f:
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
        Return skymodel names
        """
        return self.p['output_datamap']


class MergeSkymodels(Action):
    """
    Action to merge two sky models
    """
    def __init__(self, op_parset, skymodel1_datamap, skymodel2_datamap, p, prefix=None,
        direction=None, clean=True, index=None):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        skymodel1_datamap : data map
            Input data map for model 1
        skymodel2_datamap : data map
            Input data map for model 2
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
        super(MergeSkymodels, self).__init__(op_parset, 'MergeSkymodels',
            prefix=prefix, direction=direction, index=index)

        # Store input parameters
        self.skymodel1_datamap = skymodel1_datamap
        self.skymodel2_datamap = skymodel2_datamap
        self.p = p.copy()
        if self.prefix is None:
            self.prefix = 'merge_skymodel'
        self.clean = clean
        factor_working_dir = op_parset['dir_working']
        if self.direction is None:
            self.working_dir = '{0}/models/{1}/'.format(factor_working_dir,
                self.op_name)
        else:
            self.working_dir = '{0}/models/{1}/{2}/'.format(factor_working_dir,
                self.op_name, self.direction)
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        # Define script name
        self.script_file = self.parsetbasename + 'merge_skymodels.py'


    def make_datamaps(self):
        """
        Makes the required data maps
        """
        from factor.lib.datamap_lib import write_mapfile, read_mapfile

        self.p['skymodel1_datamap'] = self.skymodel1_datamap
        self.p['skymodel2_datamap'] = self.skymodel2_datamap

        model1basenames = read_mapfile(self.skymodel1_datamap)
        output_files = [os.path.splitext(bn)[0] + '_merged.skymodel' for bn in
            model1basenames]

        self.p['output_datamap'] = write_mapfile(output_files,
            self.op_name, self.name, prefix=self.prefix+'_output',
            direction=self.direction, working_dir=self.op_parset['dir_working'])


    def make_pipeline_control_parset(self):
        """
        Writes the pipeline control parset and any script files
        """
        self.p['mergescriptname'] = os.path.abspath(self.script_file)
        template = env.get_template('merge_skymodels.pipeline.parset.tpl')
        tmp = template.render(self.p)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        template = env.get_template('merge_skymodels.tpl')
        tmp = template.render(self.p)
        with open(self.script_file, 'w') as f:
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
        Return skymodel names
        """
        return self.p['output_datamap']


class FFT(Action):
    """
    Action to FFT a model into an MS
    """
    def __init__(self, op_parset, input_datamap, p, prefix=None, direction=None,
        clean=True, index=None):
        """
        Create action and run pipeline

        Parameters
        ----------
        op_parset : dict
            Parset dict of the calling operation
        input_datamap : data map
            Input data map for CASA model image
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
        super(FFT, self).__init__(op_parset, 'FFT')

        # Store input parameters
        self.input_datamap = input_datamap
        self.p = p.copy()
        if prefix is None:
            prefix = 'fft'
        self.prefix = prefix
        self.direction = direction
        self.clean = clean
        self.index = index
