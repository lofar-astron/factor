Factor: Facet Calibration for LOFAR
===================================

Factor is a tool for producing low-noise, high-resolution wide-field images from
LOFAR HBA data. Factor has been designed to use as few free parameters as
possible in order to mitigate the effects of over-fitting and thus maximize
image fidelity. Factor runs well on single machines or on compute clusters with
multiple nodes (with a shared file system). It requires only modest resources
(at least 32 GB of memory and 1 TB of disk space).


What's New in v2.0 (prerelease)
-------------------------------

* Updated to use WSClean v2.6 (including IDG support). Earlier versions are no longer supported
* Prefactor v3.0 is now required for preprocessing (earlier versions are no longer supported)
* Calibration is now done with the improved TEC solver in DPPP
* Solutions are now stored in H5parm format instead of ParmDB format; hence, LoSoTo is
now required


Installation
------------

To install Factor, follow the instructions below.

Note: Factor is currently installed on the LOFAR CEP3 cluster. Users on CEP3
should run the following command before using Factor:

    source ~rafferty/init_factor


### Dependencies

Factor requires the following:

* The LOFAR offline trunk from the LOFAR software repository (a version of the trunk after 02/11/2018 is required)
* [WSClean](http://sourceforge.net/p/wsclean/wiki/Home) (version 2.6 or later; building with IDG is recommended)
* [DP3](https://github.com/lofar-astron/DP3) (version 1.0 or later)
* [LSMTool](https://github.com/darafferty/LSMTool) (version 1.2.0 or later)
* [LoSoTo](https://github.com/revoltek/losoto) (version 2.0 or later)
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
