.. _mosaic:

Mosaic Operation
================

This operation performs the mosaicking of all facet images into a single image. The pipeline parset for this operation is ``makemosaic_pipeline.parset``.


Make mosaic
-----------

Input
	Facet images from the self-calibration (:ref:`facet_selfcal`) or the facet imaging (:ref:`facet_image`).

Output
    Single image of all facets.

Pipeline Steps
    create_images_map, create_vertices_map
        Make datamaps for input images and the facet-vertices files.

    create_compressed_mapfile_images, create_compressed_mapfile_vertices
        Create datamaps suitable for the mosaic script.

    copy_beam_info
        Copy beam information to WSClean MFS images.

    make_mosaic
        Run mosaicking script to produce final image.

Test data
    With the facet images produced by the facet self-calibration operation (:ref:`facet_selfcal`) or, if done, by the final imaging operation (:ref:`facet_image`), this step produces the mosaicked image ``RX42_SB070-079.2ch10s.wsclean2-image.make_mosaic`` in ``Test_run/results/makemosaic/field``.


Correct mosaic for primary beam
-------------------------------
Input
	Mosaic image from the previous step.

Output
    Primary-beam corrected image of all facets.

Pipeline Steps
    create_compressed_mapfile, concat
        Create datamaps suitable for the Awimager.

    make_pbimage, zero_avgpb, image2fits
        Run the Awimager to produce the ``avgpb`` image, which is then zeroed at large radii and converted to a FITS image.

    correct_mosaic
        The mosaic image is corrected with the ``avgpb`` image.

Test data
    With the mosaic image produced above, this step produces the PB-corrected image ``??`` in ``Test_run/results/makemosaic/field``.
