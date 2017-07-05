.. _runfactor:

Starting a Factor run with ``runfactor``
----------------------------------------

Factor can be run with::

    $ runfactor factor.parset

where ``factor.parset`` is the parset described in :ref:`factor_parset`. A number of options are available and are described below:

    Usage: runfactor parset

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -d                    enable dry-run mode
      -q                    enable quiet mode
      -r RESET, --reset=RESET
                            comma-separated list of directions to reset (e.g., "-r
                            facet1,facet3")
      -o OPS, --ops=OPS     comma-separated list of operations to reset for the
                            directions specified with "-r" or "--reset" (e.g., "-o
                            facetselfcal,facetsub"). By default, all operations are
                            reset. Available operations are: outlierpeel,
                            facetpeel, facetpeelimage, facetselfcal, facetsub,
                            facetimage, fieldmosaic
      -s STOP_AFTER, --stop_after=STOP_AFTER
                            Stop after processing so many facetselfcal groups.
                            (Please note the difference between directions and
                            groups!)
      -v                    enable verbose mode

Factor begins a run by checking the input measurement sets and the direction-independent instrument tables. If the instrument tables contain real/imaginary values, they are converted the phase/amplitude. The input measurement sets are chunked in time to allow more efficient processing.

Next, Factor will check for a file describing the DDE calibrators. If not found, Factor will select the DDE calibrators automatically and generate the facet regions (saved in a text file in the working directory called ``factor_directions.txt``). At this point, if ``interactive = True`` is set in the parset, Factor will pause to allow a check on the facets.

After intialization, Factor will begin self calibration and imaging (described in detail in :ref:`operations`) of the first facet or group of facets. For the facet or facets that pass self-calibration verification, Factor subtracts the improved model using the direction-dependent instrument tables. Processing then proceeds to the self calibration of the next facet or facet group, and the steps are looped until all facets have been processed.

After self calibration is finished, imaging is done for all facets. Lastly, all facet images are mosaicked together and the primary beam attenuation is corrected to produce the final image.

Factor uses the LOFAR pipeline framework to handle the actual processing. The LOFAR pipeline framework handles the distribution of jobs and keeps track of the state of a reduction. Each Factor operation is done in a separate pipeline. See :ref:`structure` for an overview of the various operations that Factor performs and their relation to one another, and see :ref:`operations` for details of each operation and their primary data products.


Checking a Factor run with ``checkfactor``
------------------------------------------

You can check the progress of a run with::

    $ checkfactor factor.parset

.. note::

    A number of ``checkfactor`` options are configurable in the Factor parset (see :ref:`parset_checkfactor_options`).


Resuming an interrupted run
---------------------------

Due to the potentially long run times and the consequent non-negligible chance
of some unforeseen failure occurring, Factor has been designed to allow easy
resumption of a reduction from a saved state and will skip over any steps that
were successfully completed previously. In this way, one can quickly resume a
reduction that was halted (either by the user or due to some problem) by simply
re-running Factor with the same parset.

For example, one can specify that only the first 5 directions be processed.
Once these directions are done, Factor will exit. One can then alter the parset
and specify that 10 directions should be done. Upon restarting, Factor will skip
over the first 5 directions and start with the 6th one (and ending with the 10th
one).

.. note::

    Upon resuming a run, Factor will pick up changes to the parset and directions file. However, changes that alter the facet layout will result in incorrect results, as Factor does not handle these properly.


Resetting a direction
---------------------

Factor allows for the processing of a direction to be reset. Resetting involves deleting the output products, resetting the state (which tracks which operations have been completed), and if necessary undoing the subtraction of the model resulting from self calibration and the imaging of the facet. Directions can be reset using the ``-r`` flag to ``runfactor``. For example, to reset all operations for direction1 and direction2::

    $ runfactor factor.parset -r direction1,direction2

Additionally, one or more specific operations can reset by including ``-o`` flag. For example, the following would reset only the facetimage operation for direction1 and direction2::

    $ runfactor factor.parset -r direction1,direction2 -o facetimage

.. note::

    If you want to reset all directions and all operations (i.e., to start the processing over from the very start), you can simply delete (or move) the ``results`` and ``state`` directories in the Factor working directory (see below), then restart Factor (without the ``-r`` or ``-o`` flags). Factor will then start
    the entire reduction again, but will skip the chunking of the input data files.
