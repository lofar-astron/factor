.. _directions_file:

The Factor directions file
==========================

A file defining the direction-dependent calibrators may be provided. The directions files is a simple text file listing various parameters of the direction-dependent calibrators (e.g., position, size, etc.). A typical directions file is shown below::

    # This is an example directions file for Factor.
    #
    # Directions should be sorted in order of reduction (e.g., bright to faint).
    #
    # Columns are defined as follows:
    #
    # name position atrous_do mscale_field_do cal_imsize solint_ph solint_amp dynamic_range region_selfcal region_facet peel_skymodel outlier_source cal_size_deg cal_flux_mjy
    #
    # Values of "empty" (for string or boolean entries) or 0 (for integer entries) indicate
    # that they should be derived internally by Factor.

    s1  14h41m01.884,+35d30m31.52 empty empty 512  1  30  LD empty empty /full/path/to/s1.skymodel True 0.01 2400
    s2  14h38m29.584,+33d57m37.82 True  True    0  1   0  LD empty empty empty False 0.02 1300
    s25 14h21m07.482,+35d35m22.87 empty empty   0  2   0  LD empty /full/path/to/s25.rgn empty False 0.1 1020
    ...

The columns are described below. The calibrators will be processed in the order in which they are given. If a directions file is provided, it must be specified with the :term:`directions_file` option in the parset. If not, Factor will identify suitable calibrators (or calibrator groups) using the values of the :term:`flux_min_for_merging_Jy`, :term:`size_max_arcmin`, :term:`separation_max_arcmin`, and :term:`flux_min_Jy` parameters in the parset.

.. note::

    You can let Factor generate the directions file itself and then edit the file to adjust the calibrators by hand. To do this, set :term:`interactive` to  ``True`` in the parset so that Factor will pause to allow editing. Once you are done, you can continue the same Factor run and it will pick up any changes made to the file.


Columns
-------

.. glossary::

    name
        The name of the calibrator or facet.

    position
        The RA and Dec (J2000) of the calibrator, written as ``RA,Dec`` in sexagesimal format (note that there should be no space after the comma). If the calibrator is a group of more than one source, this position should define the center of a calibration group.

    atrous_do
        If ``True``, the wavelet module of PyBDSM will be used during facet imaging. If empty, Factor will activate the wavelet module if it identifies a source with a diameter of 6 arcmin or larger in the facet sky model.

    mscale_field_do
        If ``True``, multiscale clean will be used during facet imaging. If empty, Factor will activate multiscale clean if it identifies a source with a diameter of 6 arcmin or larger in the facet sky model.

    cal_imsize
        The width in pixels of calibrator image. If ``0``, Factor will determine the width from the size of the calibrator or calibrator group.

    solint_ph
        The solution interval in seconds for the phase-only solve. If ``0``, Factor will set the solution interval based on the brightness of the calibrator or calibrator group.

    solint_amp
        The solution interval in seconds for the amplitude solve. If ``0``, Factor will set the solution interval based on the brightness of the calibrator or calibrator group.

    dynamic_range
        If ``HD``, amplitudes are solved for every channel. If ``LD``, amplitudes are solved in blocks defined by the :term:`TEC_block_MHz` option in the parset.

    region_selfcal
        The region to use as a clean mask during self calibration. If given, this region will be unioned with the PyBDSM-generated one.

    region_facet
        The region to use as a clean mask during facet imaging. If given, this region will be unioned with the PyBDSM-generated one.

    peel_skymodel
        The sky model to use during peeling (if the :term:`outlier_source` column is set to ``True`` or if the calibrator flux density exceeds that set with the :term:`peel_flux_Jy` option in the parset).

    outlier_source
        If ``True``, the calibrator will be peeled using the sky model given in the :term:`peel_skymodel` column and self calibration will not be done.

    cal_size_deg
        The size of the calibrator or calibrator group in degrees.

    cal_flux_mjy
        The total flux density of the calibrator or calibrator group in mJy.


