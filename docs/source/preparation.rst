.. _data_preparation:

Data preparation
================

.. note::

    Factor currently works only on HBA LOFAR data. It supports interleaved and multi-night datasets
    as well as continuous observations. Processing of international stations is not currently supported.

Factor requires that the input data be prepared using the pre-facet-calibration and initial-subtraction pipelines. These pipelines perform the calibration of the calibrator data, the removal of instrumental effects (e.g., station clock offsets), the setting of the overall amplitude scale, the calibration of the target data, and the subtraction of sources in the target field. The pipelines are available at https://github.com/lofar-astron/prefactor and must be run before Factor can be used to perform facet calibration.
