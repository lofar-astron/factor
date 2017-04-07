Factor: Facet Calibration for LOFAR
===================================

Factor is a tool for producing low-noise, high-resolution wide-field images from
LOFAR HBA data. Factor has been designed to use as few free parameters as
possible in order to mitigate the effects of over-fitting and thus maximize
image fidelity. Factor runs well on single machines or on compute clusters with
multiple nodes (with a shared file system). It requires only modest resources
(at least 32 GB of memory and 1 TB of disk space).

What's New in v1.2
------------------

* An unarchiving tool (`unarchivefactor`) has been added that can unarchive an
archive made with `archivefactor`
* An option (`update_selfcal_clean_regions`) has been added that controls
whether user-supplied clean masks are updated during selfcal
* Intersections in user-supplied clean masks are now detected and an error raised
* An archiving tool (`archivefactor`) has been added that can archive the
subtracted datasets, the sky models, the instrument tables, the selfcal plots,
and the calibrated data for one or more directions
* Polynomial sky models generated directly by WSClean during imaging are used
for prediction, resulting in improved and faster subtraction of extended
sources. Due to this change, Factor now requires WSClean v2.3 or higher


What's New in v1.1
------------------

* The subtraction of sources after self calibration has been improved
* The baseline-dependent averaging used during imaging has been reduced,
as the previous averaging caused significant smearing away from the
facet center
* Many bug fixes


Installation
------------

To install Factor, follow the instructions below.

Note: Factor is currently installed on the LOFAR CEP3 cluster. Users on CEP3
should run the following command before using Factor:

    source ~rafferty/init_factor


### Dependencies

Factor requires the following:

* The LOFAR offline trunk from the LOFAR software repository (a version of the trunk after 02/11/2016 is required)
* [WSClean](http://sourceforge.net/p/wsclean/wiki/Home) (version 2.3 or later)
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
