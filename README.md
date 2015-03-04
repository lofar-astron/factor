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

* [CASA][http://casa.nrao.org/]
* [LSMTool][http://www.numpy.org/]
* [LoSoTo][http://www.matplotlib.org/]
* [jinja2][http://jinja.pocoo.org/docs/dev/]

Usage
-----

The Factor executable can be used from the command line with a parset that defines the parameters
of the run. E.g.:

    $ factor factor.parset

The parset defines the data and working directories, etc. For details, please
see the [example parset](parsets/factor.parset).
