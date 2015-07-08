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
    TODO
