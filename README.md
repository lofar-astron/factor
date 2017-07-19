Factor: Facet Calibration for LOFAR
===================================

Factor is a tool for producing low-noise, high-resolution wide-field images from
LOFAR HBA data. Factor has been designed to use as few free parameters as
possible in order to mitigate the effects of over-fitting and thus maximize
image fidelity. Factor runs well on single machines or on compute clusters with
multiple nodes (with a shared file system). It requires only modest resources
(at least 32 GB of memory and 1 TB of disk space).


What's New Since v1.3
---------------------

* Options have been added to allow the scales used by WSClean during
clean to be specified (`selfcal_multiscale_scales_pixel` and
`facet_multiscale_scales_pixel` under `[imaging]`)

* The `atrous_do` column in the directions file has been renamed to
`mscale_selfcal_do`, as this name better describes its purpose

* Improved handling of flagged solutions during smoothing of the amplitude
solutions


What's New in v1.3
------------------

* Updated to use WSClean 2.4. Earlier versions are no longer supported
* Preaveraging is now done in frequency as well as in time. This
preaveraging generally improves the S/N of the CS slow-gain solutions for
fainter sources by a factor of ~ 2
* An option (`min_fraction_per_band`) has been added that sets the
minimum allowed unflagged fraction per band
* WSClean's automasking feature is now used during imaging. The old
image-mask-image sequence is no longer used during self calibration, but can
still be used during the final, full-bandwidth facet imaging if
`automask_facet_image = False` under the `[imaging]` section of the parset. The
option `update_selfcal_clean_regions` has been removed, as it
no longer applies (PyBDSF masking has been removed from selfcal)
* The combination of flagging ranges specified by the `flag_abstime`,
`flag_baseline`, and `flag_freqrange` options can now be set with the
`flag_expr` option


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


Installation
------------

To install Factor, follow the instructions below.

Note: Factor is currently installed on the LOFAR CEP3 cluster. Users on CEP3
should run the following command before using Factor:

    source ~rafferty/init_factor


### Dependencies

Factor requires the following:

* The LOFAR offline trunk from the LOFAR software repository (a version of the trunk after 02/11/2016 is required)
* [WSClean](http://sourceforge.net/p/wsclean/wiki/Home) (version 2.4 or later)
* [LSMTool](https://github.com/darafferty/LSMTool) (version 1.2.0 or later)
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
