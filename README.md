FACTOR: Facet Calibration for LOFAR
===================================

Installation
------------

FACTOR is currently installed on the LOFAR CEP3 cluster. Users on CEP3
should run the following command before using FACTOR:

    source ~rafferty/init_factor

If you want to install FACTOR yourself, follow the instructions below.

### Dependencies

FACTOR requires the following:

* The LOFAR offline trunk from the LOFAR software repository (version 2.15 or later)
* [CASA](http://casa.nrao.edu) (versions *before* 4.3 or version 4.5.0 or later; note that CASA requires Java for the buildmytasks script that is needed by FACTOR, so you will need Java as well)
* [WSClean](http://sourceforge.net/p/wsclean/wiki/Home) (version 1.10 or later)
* [LSMTool](https://github.com/darafferty/LSMTool) (version 1.1 or later)
* [jinja2](http://jinja.pocoo.org/docs/dev)
* [Shapely](https://github.com/Toblerity/Shapely)
* [APLpy](http://aplpy.github.io) (version 1.0 or later)

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

FACTOR may be thought of as a tool to set up and run LOFAR pipeline scripts that
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
* `parset.py`: handles the reading of the FACTOR parset
* `process.py`: performs the actual processing (running of pipelines, etc.)

Usage
-----

The FACTOR executable (named `runfactor`) can be used from the command line with
a parset that defines the parameters of the run. E.g.:

    $ runfactor factor.parset

You can check the progress of a run with `checkfactor`:

    $ checkfactor factor.parset

The parset defines the data and working directories, various options, etc.
FACTOR handles all the initialization and sets up the directories, pipeline
parsets, etc. FACTOR will also run the pipelines by default. However, if you
want to run the pipelines yourself, you can set the dry_run (`-d`) argument
(e.g., `runfactor -d factor.parset`). FACTOR will then perform all the setup but
will not run the pipelines. Note that this option is most useful if a
directions file is given, as without this file FACTOR cannot create pipelines
for self calibration.

For details on the usage, please see the [full documentation](http://www.astron.nl/citt/facet-doc/)
and the examples in the examples directory.
