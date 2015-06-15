Factor: Facet Calibration for LOFAR
===================================

Installation
------------

Factor is currently installed on the LOFAR CEP3 cluster. Users on CEP3
should run the following commands before using Factor:

    use LofIm
    source ~rafferty/init_factor

### Requirements

Beyond the standard LOFAR software packages, Factor requires the following:

* [CASA][http://casa.nrao.org]
* [WSClean][http://sourceforge.net/p/wsclean/wiki/Home]
* [LSMTool][https://github.com/darafferty/LSMTool]
* [LoSoTo][https://github.com/revotek/losoto]
* [jinja2][http://jinja.pocoo.org/docs/dev]
* [Shapely][https://github.com/Toblerity/Shapely]

Usage
-----

The Factor executable can be used from the command line with a parset that
defines the parameters of the run. E.g.:

    $ factor factor.parset

The parset defines the data and working directories, all options, etc. For
details, please see the full documentation and the examples in the examples
directory.

Input data
----------

Currently, only HBA data are supported. The input MS files should be
concatenated bands of 10 subbands. They should have clock offsets removed, the
beam at the phase center applied (to DATA and CORRECTED\_DATA columns), and be
calibrated to 25 arcsec resolution (in the CORRECTED\_DATA column). All MS files
should be placed together in a single directory.
