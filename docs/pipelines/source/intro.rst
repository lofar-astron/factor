.. _intro:

Introduction
============

The Factor (Facet Calibration for LOFAR) software package has been developed to perform facet calibration in a user-friendly way for LOFAR HBA data. Running Factor on a dataset generally involves the following steps:

Preparing the data (see :ref:`data_preparation`)
    Factor requires that the input data be prepared in a specific way. This step involves running the pre-Factor pipelines.

Making a Factor parset (see :ref:`factor_parset`)
    This step involves making the parset that defines the options of the Factor run.

Making a directions file (see :ref:`directions_file`)
    (Optional) This step involves making a file that defines the directions toward which direction-dependent solutions will be derived.

Running Factor (see :ref:`running_factor`)
    This step involves running Factor and checking its progress and output.
