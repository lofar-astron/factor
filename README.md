Factor: Facet Calibration for LOFAR
===================================

Installation
------------

Factor is currently installed on the LOFAR CEP3 cluster. Users on CEP3
should run the following commands before using Factor:

    use LofIm
    source ~rafferty/init_factor

### Dependencies

Beyond the standard LOFAR software packages, Factor requires the following:

* [CASA](http://casa.nrao.org) (versions *before* 4.3)
* [WSClean](http://sourceforge.net/p/wsclean/wiki/Home)
* [LSMTool](https://github.com/darafferty/LSMTool)
* [jinja2](http://jinja.pocoo.org/docs/dev)
* [Shapely](https://github.com/Toblerity/Shapely)

Code Structure
--------------
The top-level directories are:

* `bin`: contains the `runfactor` executable
* `docs`: contains the documentation
* `examples`: contains examples parset and directions file
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

* `cluster.py`: handles the compute cluster setup
* `directions.py`: handles the setup of directions (facets)
* `parset.py`: handles the reading of the Factor parset
* `process.py`: preforms the actual processing (running of pipelines, etc.)

Usage
-----

The Factor executable can be used from the command line with a parset that
defines the parameters of the run. E.g.:

    $ runfactor factor.parset

Alternatively, you can run factor from Python scripts with:

    from factor import process
    process.run('factor.parset')

The parset defines the data and working directories, various options, etc.
Factor handles all the initialization and sets up the directories, pipeline
parsets, etc. Factor will also run the pipelines by default. However, if you
want to run the pipelines yourself, you can set the dry_run (`-d`) argument
(e.g., `runfactor -d factor.parset`). Factor will then perform all the setup but
will not run the pipelines. Note that this option is most useful if a
directions file is given, as without this file Factor cannot create pipelines
for self calibration, etc.

For details on the usage, please see the full documentation and the examples in
the examples directory.
