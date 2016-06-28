Factor: Facet Calibration for LOFAR
===================================

Installation
------------

Factor is currently installed on the LOFAR CEP3 cluster. Users on CEP3
should run the following command before using Factor:

    source ~rafferty/init_factor

If you want to install Factor yourself, follow the instructions below.

### Dependencies

Factor requires the following:

* The LOFAR offline trunk from the LOFAR software repository (version 2.16.4 or later)
* [CASA](http://casa.nrao.edu) (version 4.2 or version 4.5 -- versions 4.3, 4.4, and 4.6 do not work; note that CASA requires Java for the buildmytasks script that is needed by Factor, so you will need Java as well)
* [WSClean](http://sourceforge.net/p/wsclean/wiki/Home) (version 1.10 or later)
* [LSMTool](https://github.com/darafferty/LSMTool) (version 1.1 or later)
* [jinja2](http://jinja.pocoo.org/docs/dev)
* [Shapely](https://github.com/Toblerity/Shapely)
* [APLpy](http://aplpy.github.io) (version 1.0 or later)
* [pyds9](https://github.com/ericmandel/pyds9) (optional, to allow checkfactor to interface with ds9)

### Downloading and Installing

Get the latest developer version by cloning the git repository:

    git clone https://github.com/lofar-astron/factor.git

Then install with:

    cd factor
    python setup.py install

Code Structure
--------------
The top-level directories are:

* `bin`: contains the `runfactor`  and `checkfactor` executables
* `docs`: contains the documentation
* `examples`: contains example parset and directions file
* `factor`: main code tree

Factor may be thought of as a tool to set up and run LOFAR pipeline scripts that
perform facet calibration. As such, the main code tree is comprised of the
following directories:

* `lib`: contains general classes used to define operations, etc.
* `operations`: contains the operations to perform selfcal, etc.
* `parsets`: contains parsets used by executables called by the pipeline, such
as DPPP, etc.
* `pipeline`: contains all the pipeline files, such as pipeline parsets,
configuration files, and task definitions
* `scripts`: contains custom scripts called by the pipeline
* `skymodels`: contains generic sky models needed by the pipeline

In addition, the `factor` directory contains the following files:

* `check_progress.py`: runs the interactive progress tool
* `cluster.py`: handles the compute cluster setup
* `directions.py`: handles the setup of directions (facets)
* `parset.py`: handles the reading of the Factor parset
* `process.py`: performs the actual processing (running of pipelines, etc.)

Usage
-----

The Factor executable (named `runfactor`) can be used from the command line with
a parset that defines the parameters of the run. E.g.:

    $ runfactor factor.parset

You can check the progress of a run with `checkfactor`:

    $ checkfactor factor.parset

The parset defines the data and working directories, various options, etc.
Factor handles all the initialization and sets up the directories, pipeline
parsets, etc.

For details on the usage, please see the [full documentation](http://www.astron.nl/citt/facet-doc/)
and the examples in the examples directory.
