Initial Subtract Operation
==========================

This operation images each band at high and low resolution to make and subtract
sky models.

.. _high_res_images:

Make high-res images
--------------------

Input
	MS files concatenated into bands of 10 subbands (2 MHz) each.

	The data should have clock offsets removed, the beam at the phase center
	applied (to ``DATA`` and ``CORRECTED_DATA`` columns), and be calibrated to
	25 arcsec resolution (in the ``CORRECTED_DATA`` column).

Output
    For each band, a wide-field (~ 7 degree radius) image is made at
    approximately 25 arcsec resolution.

    .. note::

        The image should show the
        higher-resolution structure. There may be strong artifacts around bright
        sources.

    A small region (~ 1 degree square) of an example image is shown
    below.

    .. figure:: high_res_image.png
       :scale: 40 %
       :figwidth: 75 %
       :align: center
       :alt: example image

Test data
    TODO

Make high-res sky models
------------------------

Input
    Model images created during high-res imaging (see :ref:`high_res_images`), one for each band.

Output
    Sky models in ``makesourcedb`` format for each band

Test data
    TODO
