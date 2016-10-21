Factor: Facet Calibration for LOFAR
===================================

Factor is a tool for producing low-noise, high-resolution wide-field images from LOFAR HBA data. Factor has been designed to use as few free parameters as possible in order to mitigate the effects of over-fitting and thus to maximize image fidelity. Factor runs well on single machines or on compute clusters with multiple nodes. It requires only modest resources (at least 32 GB of memory and 1 TB of disk space).

What's New
----------

Note: Factor v1.0 is not compatible with previous versions. If you started a reduction
with a previous version, you must restart from scratch to use v1.0; the version pervious to v1.0 is still
available on the "pre1.0" branch.

Since version pre1.0:

* The BBS calibration software has been replaced by GainCal, which is both much faster and more stable
* The CASA imager has been replaced by WSClean, which is generally much faster and integrates better with the LOFAR pipeline framework that Factor uses
* A new option has been added to keep the primary data files used in self calibration in memory, speeding up self calibration dramatically on some systems (set "dir\_local\_selfcal = /dev/shm" under the [cluster] section of the parset)
* The Dysco storage manager can now be used to compress visibilities and weights, reducing file sizes (and IO) by a factor of ~ 2.5
* The reimage\_selfcaled option has been removed, as facets are always reimaged now. A image\_target\_only option has been added if you wish to reimage only the target
* The skip\_facet\_imaging option has been removed, as facets are always imaged now to improve the subtraction (using up to 6 bands)


Installation
------------

To install Factor, follow the instructions below.

Note: Factor is currently installed on the LOFAR CEP3 cluster. Users on CEP3
should run the following command before using Factor:

    source ~rafferty/init_factor


### Dependencies

Factor requires the following:

* The LOFAR offline trunk from the LOFAR software repository (release 2.18 or later; to use compression, a version of the trunk after 14/10/2016 is required)
* [WSClean](http://sourceforge.net/p/wsclean/wiki/Home) (version 1.12 or later)
* [LSMTool](https://github.com/darafferty/LSMTool) (version 1.1 or later)
* [jinja2](http://jinja.pocoo.org/docs/dev)
* [Shapely](https://github.com/Toblerity/Shapely)
* [APLpy](http://aplpy.github.io) (version 1.0 or later)
* [pyds9](https://github.com/ericmandel/pyds9) (optional, to allow checkfactor to interface with ds9)
* [Dysco](https://github.com/aroffringa/dysco) (optional, to allow compression of the visibilities and weights)

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
