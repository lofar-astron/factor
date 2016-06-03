.. _capabilities:

Capabilities of Factor
======================

Factor corrects for direction-dependent effects in HBA LOFAR data, including ionospheric effects and beam-model errors. These corrections are essential to obtaining instrumental-noise limited (~ 0.1 mJy/beam for an 8-hour observation), high-resolution (~ 5 arcsec FWHM) images. Factor works by dividing up the field into many facets and solving for the direction-dependent corrections in each facet. Importantly, Factor was designed to minimize the number of free parameters needed to parameterize these corrections to avoid overfitting. This minimization is critical to producing high-fidelity images.

Factor supports interleaved and multi-night datasets as well as continuous observations. Factor was designed to work on compute clusters and allows for the distribution of jobs over multiple nodes of a cluster and for the processing of facets in parallel.

