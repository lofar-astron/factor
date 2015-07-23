.. _pre_facet:

Pre-facet Calibration Pipeline
==============================

This operation performs a direction independent calibration of the target data following the procedure outlined in van Weeren et al (in prep) and Williams et al (in prep) and users are strongly encouraged to read those papers before using this pipeline. First the amplitude and clock calibration solutions are determined from a bright calibrator source and then these solutions are applied to the target field which is finally phase calibrated using a sky model. In the following subsections each part of this direction independent calibration pipeline parset is detailed.

.. warning::
    The treatment of the LOFAR beam in this pipeline is
    correct if imaging the final products in awimager to produce images
    with the correct application of the LOFAR beam but without direction
    dependent calibration. However, if using the facet calibration
    pipeline there are alterations required to how the LOFAR beam is
    applied. These alterations are relevant to Sections
    :ref:`section_calib`, :ref:`sec:parmmap_calibtarget` and
    :ref:`sec:gsmcalibtarget` and are described in Section
    :ref:`pre-facet-calibration-alterations`.


.. _section_timmap:

timmap
------

In this stage we create a list of the calibrator files that will be processed in later steps.

timmap.control.opts.folder = {{user_caldata_dir}}
    The full path to a directory that contains only the calibrator measurement sets which have uncalibrated data in the DATA column.


.. _section_dppprepcal:

dpppprepcal
-----------

In this step we aim to remove bad data, international stations and also average the calibrator data (if necessary). The output of Section :ref:`section_timmap` is the input for this NDPPP operation. NDPPP parameters that may require alteration are detailed below.

dpppprepcal.argument.avg.freqstep = 2
    Average the calibrator dataset in frequency. In this example the initial dataset had an averaging of 8ch/sb so this averaging produces a dataset with an averaging of 4ch/sb.
dpppprepcal.argument.avg.timestep  = 2
    Average the calibrator dataset in time. In this example the initial dataset had an averaging of 2seconds so this averaging produces a dataset with an averaging of 4seconds.
dpppprepcal.argument.flag.baseline = [ CS013HBA* ]
    Flag stations where the data are bad (the dipoles on CS013* are rotated so the data should be flagged).
dpppprepcal.argument.flagint.baseline = UK*&&*; SE*&&* ; FR*&&*; DE*&&*
    Flag international stations.
dpppprepcal.argument.filterint.baseline = !UK*&&* ; !SE*&&* ; !FR*&&* ; !DE*&&*
    Remove international stations from the measurement set to reduce its size.
dpppprepcal.control.opts.max_per_node = 6
    The number NDPPP runs to do simultaneously

The output files from this NDPPP are located in  ``working_dir`` (which is specified in ``pipeline.cfg`` setup file) followed by the pipeline parset name.

.. _section_calibcal:

calibcal
--------
The flagged and averaged measurement sets output from Section :ref:`section_dppprepcal` are now calibrated with BBS using a model of the calibrator.

calibcal.argument.parset = {{ user_parset_path }}calibcal.parset
    The BBS parset that will be used to calibrate the calibrator.
calibcal.argument.catalog = {{ user_skymodel_path }}3C196-pandey.skymodel
    The appropriate sky model for the calibrator.
calibcal.control.max_per_node = 6
    The number of BBS runs to do simultaneously.

The parset used for this BBS calibration is given below. This parset determines and applies one solution per channel for the XX and YY gains as well as the rotation angle to take account of differential Faraday Rotation::

    Strategy.InputColumn = DATA
    Strategy.TimeRange = []
    Strategy.Baselines = *&
    Strategy.ChunkSize = 1000
    Strategy.UseSolver = F
    Strategy.Steps = [solve, correct]

    Step.solve.Operation = SOLVE
    Step.solve.Model.Sources = []
    Step.solve.Model.Beam.Enable = True
    Step.solve.Model.Cache.Enable = T
    Step.solve.Model.Phasors.Enable = F
    Step.solve.Model.Bandpass.Enable = F
    Step.solve.Model.Gain.Enable = T
    Step.solve.Model.DirectionalGain.Enable = F
    Step.solve.Model.CommonRotation.Enable = T
    Step.solve.Solve.Parms = ["Gain:0:0:*","Gain:1:1:*","CommonRotationAngle:*"]
    Step.solve.Solve.UVRange = [160,50000]
    Step.solve.Solve.ExclParms = []
    Step.solve.Solve.CellSize.Freq = 1
    Step.solve.Solve.CellSize.Time = 1
    Step.solve.Solve.CellChunkSize = 100
    Step.solve.Solve.PropagateSolutions = T
    Step.solve.Solve.Options.MaxIter = 100
    Step.solve.Solve.Options.EpsValue = 1e-9
    Step.solve.Solve.Options.EpsDerivative = 1e-9
    Step.solve.Solve.Options.ColFactor = 1e-9
    Step.solve.Solve.Options.LMFactor = 1.0
    Step.solve.Solve.Options.BalancedEqs = F
    Step.solve.Solve.Options.UseSVD = T

    Step.correct.Operation = CORRECT
    Step.correct.Model.Sources = []
    Step.correct.Model.Gain.Enable = T
    Step.correct.Model.Beam.Enable = T
    Step.correct.Model.CommonRotation.Enable = T
    Step.correct.Output.Column = CORRECTED_DATA

The outputs of this step are instrument tables within each of the measurement sets that BBS has operated on. Diagnostic plots showing the quality of the solutions are produced later in the pipeline (see Sections  :ref:`section_fitclock`, :ref:`section_ampl` and :ref:`section_plots`).


ginst, gvdsfile, parmcoll and h5imp
-----------------------------------

Here we create a H5parm file that contains the instrument tables from all the calibrator files. This is procedure is done in several stages.

* The ginst function creates vds files for each instrument table and saves these into your ``pre-facet-calibration.output.job_dir+/vds`` directory.
* The gvdsfile function collects all the vds files to create a gvds file.
* The parmcoll function uses the LoSoTo script ``parmdb_collector.py`` to gather together all of your calibrator instrument tables into the ``workingdir_dir``.
* The h5imp function uses the LoSoTo script ``H5parm_importer.py`` to create a H5parm file that contains information from all of the instrument tables for the calibrator.

For these four functions only the paths to the LoSoTo scripts and the storage directory for your vds files should be altered.

gvdsfile.control.opts.folder
    Output folder

parmcoll.control.opts.executable
    Path to ``parmdb_collector.py``

h5imp.control.opts.executable
    Path to ``H5parm_importer.py``


.. _section_fitclock:

fitclock
--------

The H5parm file that contains the instrument tables from all of the calibrators is used to fit the clock and the TEC (see van Weeren et al., and Williams et el., for details). Only the path to the fitting script should to be altered.

h5imp.control.opts.executable={{user_script_dir}}fit_clocktec_initialguess_losoto.py
    The path to the LoSoTo script

This function outputs four numpy arrays to your ``workingdir_dir`` which contain the clock and dTEC values derived from the fitting.  Two of these arrays are for the clock values  (``fitted_data_dclock_somename_1st.npy`` and ``fitted_data_dclock_somename_1st.sm.npy``) and the other two (``fitted_data_dclock_somename_1st.sm.npy`` and ``fitted_data_dclock_somename_1st.sm.npy``) are the dTEC values. Diagnostic plots from the clock and dTEC fits will be created in Section :ref:`section_plots`.


.. _section_ampl:

ampl
----

The H5parm file that contains the instrument tables from all of the calibrators is used to determine a single XX and YY amplitude calibration value for each sub band. Only the path to the script should be altered.

ampl.control.opts.executable={{user_script_dir}}amplitudes_losoto.py
    The path to the LoSoTo script

This function outputs a numpy array containing the amplitude values as a function of frequency (``freqs_for_amplitude_array.npy``) and diagnostic plots which should all be inspected. Descriptions and examples of these plots are shown in Figure :num:`figure-amp`  and :num:`figure-amp2`.

.. _figure-amp:

.. figure:: CS031HBA0_ampmat_both.pdf
   :figwidth: 90 %

   An example of an ``*_ampmat.pdf`` plot is shown on the left and the corresponding ``*_ampmat_smooth.pdf`` plot is shown on the right. These plots show the amplitude in colour scale as a function of time and frequency. The user should search for unexpected ripples, troughs, peaks or outlier antennas which can then be identified as bad data and flagged before the pipeline and calibration is restarted. These plots are shown for each antenna in the ``matrix_xx.png`` and ``matrix_yy.png`` images.


.. _figure-amp2:

.. figure:: CS031HBA0_profileXX_YY.pdf
   :figwidth: 90 %

   An example of an ``*_profileXX.pdf`` plot is shown on the left an the corresponding ``*_profileYY.pdf`` plot is shown on the right. These plots show the derived  average calibration for the antenna as a function of frequency. The user should search for unexpected ripples, troughs, peaks or outlier antennas which can then be identified as bad data and flagged before the pipeline and calibration is restarted.


.. _section_plots:

plots
-----

Now that the amplitude, clock and dTEC solutions have all been determined from the calibrator instrument tables we run a simple script to produce some further diagnostic plots. Only the path to the script should be altered.

ampl.control.opts.executable={{user_script_dir}}examine_npys.py
    The path to the plotting script

This function outputs the diagnostic plots which are described in Figures :num:`figure-plots1`, :num:`figure-plots2` and :num:`figure-plots3`.

.. _figure-plots1:

.. figure:: dtec_allsols.png
   :figwidth: 90 %

   An example ``dtec_allsols.png`` plot which shows the derived dTEC values for each antenna as a function of time. These values are not applied to the target data because the TEC for the calibrator is different to the target. However, it is information to see the differential TEC values and their variation in time as this gives a measure of the behaviour of the ionosphere during the calibrator observation (see van Weeren et al in prep).

.. _figure-plots2:

.. figure:: dclock_allsols.png
   :figwidth: 90 %

   An example ``dClock_allsols.png`` plot which shows the derived dTEC values for each antenna as a function of time. This the difference in clock values between each antenna and the reference antenna (CS001). The core stations all have low differences in clock values whereas the ones for the remote stations are up to 100ns. The clock values should be approximately constant across the calibrator observation (see van Weren et al in prep).

.. _figure-plots3:

.. figure:: amp_allsols.png
   :figwidth: 90 %

   An example ``amp_allsols.png`` plot which shows the derived amplitude solutions for each antenna as a function of time. Each antenna should have similar amplitude calibration solutions and outlier antennas or bad sub bands can be easily spotted on this plot. Any bad data that is recognise can flagged before the pipeline and calibration is restarted.


.. _section_phase:

concatmapcal and phase
----------------------

The final stage of calibration of the calibrator is to determine a median phase offset between the XX and YY per antenna and produce a further diagnostic plot. Only the path to the script should be altered.

phase.control.opts.executable={{user_script_dir}}find_cal_global_phaseoffset.py
    The path to the script

The outputs of this script are the numpy array ``somename_phase_array.npy`` and corresponding plot ``phase_xx_yy_offset.png``. The plot is described in Figure :num:`figure-phase1` and must be inspected by the user. The file ``freqs_for_phase_array.npy`` is also output and this contains the frequencies for each sub band in your calibrator measurement sets but do not need to be inspected by the user.


.. _figure-phase1:

.. figure:: phase_xx_yy_offset.png
   :figwidth: 90 %

   An example ``phase_xx_yy_offset.png`` plot which shows the derived XX and YY phase offsets (phase on the y-axis and sub band on the x-axis). The blue shows the unsmoothed offsets and may have small peaks or trough where the values were not well determined. The green shows the smoothed values and no peaks or troughs are expected here. Any bad data that is recognise can flagged before the pipeline and calibration is restarted.


.. _section_calib:

calib
-----

.. warning::
    Please see Section :ref:`pre-facet-calibration-alterations` if using this pipeline to prepare for facet calibration.

Calibration of the calibrator is now complete and the calibration solutions are stored in numpy arrays that were created in previous steps. To apply these numpy arrays to the target data we first need to create a template parmdb instrument table which has the appropriate structure. This needs to be done for only one of your target measurement sets and only the path to your chosen target measurement set needs to altered in the parset.

calib.control.opts.arguments
    The arguments given to BBS to create an instrument table called template for the file ``L258233_SB355_uv.dppp.MS`` using the parset ``calibtarget_beam.parset`` and the sky model ``template_parmdb.model``. ``L258233_SB355_uv.dppp.MS`` can be any of your target measurement sets.

The ``template_parmdb.model`` is any sky model which contains one or more sources, it does not matter what sources are in the sky model as its purpose is just to create a template rather than to derive calibration parameters. The sky model should be kept as simple as possible to create the template parmdb at minimal computational expense. The  parset used for BBS calibration is given below. This parset creates a template parmdb that contains entries for the XX and YY gains as well as the clock. The parmdb contains one entry for each sub band and for each time chunk and the MaxIter is deliberately set to 1 for faster processing::

    Strategy.InputColumn = DATA
    Strategy.ChunkSize   = 200
    Strategy.UseSolver   = F
    Strategy.Steps       = [solve]

    Step.solve.Model.Sources                = [] # all in skymodel
    Step.solve.Model.Cache.Enable           = T
    Step.solve.Model.Phasors.Enable         = F
    Step.solve.Model.DirectionalGain.Enable = F
    Step.solve.Model.Gain.Enable            = T
    Step.solve.Model.Clock.Enable           = T
    Step.solve.Model.TEC.Enable             = F
    Step.solve.Operation                    = SOLVE
    Step.solve.Solve.Parms                  = ["Gain:0:0:*","Gain:1:1:*","Clock:*"]
    Step.solve.Solve.CellSize.Freq          = 0
    Step.solve.Solve.CellSize.Time          = 1
    Step.solve.Solve.CellChunkSize          = 100
    Step.solve.Solve.PropagateSolutions     = T
    Step.solve.Solve.Options.MaxIter        = 1
    Step.solve.Solve.Options.EpsValue       = 1e-9
    Step.solve.Solve.Options.EpsDerivative  = 1e-9
    Step.solve.Solve.Options.ColFactor      = 1e-9
    Step.solve.Solve.Options.LMFactor       = 1.0
    Step.solve.Solve.Options.BalancedEqs    = F
    Step.solve.Solve.Options.UseSVD         = T
    Step.solve.Solve.Mode                   = COMPLEX
    Step.solve.Model.Beam.Enable            = T
    Step.solve.Model.Beam.UseChannelFreq    = T


.. _section_timtargetmap:

timtargetmap
------------

In this stage we create a list of the target files that will be processed in later steps.

timtargetmap.control.folder
    The full path to a directory that contains only the target measurement sets which have uncalibrated data in the DATA column.


.. _section_dpppreptar:

dppppreptar
-----------

In this step we aim to remove bad data, international stations and also average the target data (if necessary). The output of Section :ref:`section_timtargetmap` is the input for this NDPPP operation. NDPPP parameters that may require alteration are detailed below.

dppppreptar.argument.avg.freqstep = 2
    Average the target dataset in frequency. In this example the initial dataset had an averaging of 8ch/sb so this averaging produces a dataset with an averaging of 4ch/sb.
dppppreptar.argument.avg.timestep  = 2
    Average the target dataset in time. In this example the initial dataset had an averaging of 2seconds so this averaging produces a dataset with an averaging of 4seconds.
dppppreptar.argument.flag.baseline = [ CS013HBA* ]
    Flag stations where the data are bad (the dipoles on CS013* are rotated so the data should be flagged).
dppppreptar.argument.flagint.baseline = UK*&&*; SE*&&* ; FR*&&*; DE*&&*
    Flag international stations.
dppppreptar.argument.filterint.baseline = !UK*&&* ; !SE*&&* ; !FR*&&* ; !DE*&&*
    Remove international stations from the measurement set to reduce its size.
dppppreptar.control.opts.max_per_node = 6
    The number of NDPPP runs to do simultaneously

The output files from this NDPPP are located in  ``working_dir`` (which is specified in ``pipeline.cfg`` setup file) followed by the pipeline parset name.


ateamtarget and ateamcliptar
----------------------------

Whether or not you should do this step depends on your demixing setup in your observation. The contribution of 'ateam' sources should be removed from your data to obtain lower noise levels and higher fidelity images. Here we assume that data has had no demixing and we use the ateamclipper.py to minimise the effects of the 'ateam' sources. Only the path of the BBS parset and the 'ateam' sky model needs altering.

ateamtarget.control.opts.arguments
    The arguments given to BBS to predict the response from the 'ateam' sources. The paths to the BBS parset and the sky model must be updated.
ateamtarget.control.opts.max_per_node=6
    The number of BBS runs to do simultaneously.

The BBS ``{{user_parset_dir}}ateamclip.parset`` that is used to predict the 'ateam' sources::

    Step.predict4.Model.Sources         = [VirA_4_patch,CygAGG,CasA_4_patch,TauAGG]
    Step.predict4.Model.Cache.Enable    = T
    Step.predict4.Model.Gain.Enable     = F
    Step.predict4.Operation             = PREDICT
    Step.predict4.Output.Column         = MODEL_DATA
    Step.predict4.Model.Beam.Enable     = T
    Step.predict4.Model.Beam.UseChannelFreq = T

The BBS run fills in the ``MODEL_DATA`` column with a prediction of the 'ateam' sources behaviour. Once the BBS run is completed the ``Ateamclipper.py`` script is used to remove severely contaminated data.

ateamcliptar.control.opts.executable={{user_script_dir}}ateamclipper.py
    The path to ateamclipper.py.
ateamcliptar.control.opts.max_per_node=6
    The number of ateamclipper.py runs to do simultaneously.


.. _section_trans:

trans
-----

Now that the target data have been flagged and the contribution from 'ateam' sources has been minimised we can transfer the calibration solutions from the calibrator to correct for clock offsets, phase offsets and amplitude. To transfer these values we use the template instrument table that was created in Section :ref:`section_calib` and copy this to each measurement set of the target and fill it up with the appropriate values from the calibration numpy arrays created in Sections :ref:`section_fitclock`, :ref:`section_ampl` and :ref:`section_phase`. Only the paths to the datasets and scripts needs to be altered.

trans.control.opts.executable={{user_scrip_dir}}transfer_amplitudes+clock+offset.py
    The path to the script which transfers solutions from numpy arrays to instrument tables
trans.control.opts.arguments
    The arguments required for the script. Here ``instrument_amp_clock`` is the name of the instrument table that will be created for each measurement set of you target. ``{{working_dir}}{{parset_name}}`` is the name of the output directory of the pipeline and ``{{user_data_path}L258233_SB355_uv.dppp.MS/template`` is the template instrument table that was created in Section :ref:`section_calib`.

The outputs of this step are ``instrument_amp_clock_offset`` tables in each of your target measurement sets in your ``{{working_dir}}``.



.. _sec:parmmap_calibtarget:

parmmap and calibtarget
-----------------------

.. warning::
    Please see Section :ref:`pre-facet-calibration-alterations` if using this pipeline to prepare for facet calibration.

In Section :ref:`section_trans` we created the appropriate instrument tables to correct the clock, phase offset and amplitude of our target data and in this section we apply those instrument tables. The function ``parmmap`` is used to create a new list of files that contains all of the ``instrument_amp_clock_offset`` tables that were made in Section :ref:`section_trans` and none of the inputs to this function need altering. In the ``calibtarget`` we apply the instrument tables and the user must ensure that the paths are correct.

calibtarget.control.opts.arguments=[-v,--parmdb,ampinstrument,targetms, {{user_parset_directory}applyparmdb.parset]
    The arguments given to BBS to apply the appropriate instrument table. The path to the ``applyparmdb.parset`` must be given by the user.
calibtarget.control.opts.max_per_node=6
    The number of BBS runs to do simultaneously.

The BBS parset which applies the existing instrument tables is given below::

    Strategy.InputColumn = DATA
    Strategy.ChunkSize = 1000
    Strategy.UseSolver = F
    Strategy.Steps = [ correct]

    Step.correct.Operation = CORRECT
    Step.correct.Model.Sources = []
    Step.correct.Model.Cache.Enable  = T
    Step.correct.Model.Clock.Enable = T
    Step.correct.Model.Gain.Enable = T
    Step.correct.Model.CommonRotation.Enable = F
    Step.correct.Model.Beam.Enable = F
    Step.correct.Model.Beam.UseChannelFreq = T
    Step.correct.Output.Column = CORRECTED_DATA


dpppaverage
-----------

Once the calibration has been applied and the ateam sources have been demixed we can flag the data again and also average a little further if permitted by time average smearing and frequency averaging smearing limits (see e.g. Bridle & Schwab (1989) for approximate formulas to calculate the effects of smearing). Only the averaging and the number of simultaneous runs need to be altered by the user.

dpppaverage.control.opts.max_per_node = 6
    The number of NDPPP runs to do simultaneously.
dpppaverage.parsetarg.avg.freqstep = 2
    Desired frequency averaging.
dpppaverage.parsetarg.avg.timestep = 2
    Desired time averaging.


conatmaptar, createmap2 and dpppconcat
--------------------------------------

Once the data has been averaged we combine into datasets that contain more than one sub band before phase calibration.. To do this we use the ``concatmaptar`` and ``createmap2`` functions which combined create a map file that links a given number of sub bands into one output file. We then use ``dpppconcat`` to combine the sub bands. In this section the only alteration required is the number of sub bands that the user combines.

createmap2.control.opts.listsize=12
    The number of sub bands to combine into one measurement set prior to phase calibration. The appropriate number to combine depends a little on the observing conditions (see van Weeren et al in prep) and here we use 12.


flagrfi
-------

After the data has been combined we search for RFI again in an attempt to pick up low level RFI that was previously missed when just examining individual sub bands. In this function nothing needs to be altered by the user.


.. _sec:gsmcalibtarget:

gsmcalibtarget
--------------

.. warning::
    Please see Section :ref:`pre-facet-calibration-alterations` if using this pipeline to prepare for facet calibration.

The final step in the direction independent calibration of the target field is to calibrate the data off an existing sky model. Here we use a model generated by the gsm.py but any appropriate model can be used. The user must ensure the path to the sky model is correct and also that the number of simultaneous jobs is suitable.

gsmcalibtarget.control.opts.arguments
    The arguments given to BBS to phase calibrate the data. The path to past and the model (here ``P21_hetdex_5deg.txt``) must be given by the user.
gsmcalibtarget.control.opts.max_per_node = 6
    The number of BBs jobs to run simultaneously.

The BBS parset used to calibrate the phase of the target data is given below. A single solution is found for all frequencies in the merged sub bands and for each 4 time samples. It is best to keep the averaging of the data low so that a solution is found on approximately a 20-30 second timescale to ensure there is sufficient signal but that the ionosphere does not vary significantly within this timescale. Example phase solutions are shown in Figure :num:`figure-phase`. Time slots or antennas with noticeable phase errors (where the phase behaviour changes rapidly and looses its structure or antennas with consistently incoherent phases) can be flagged after the phase calibration is complete. The parset is::

    Strategy.InputColumn = DATA
    Strategy.ChunkSize = 200
    Strategy.UseSolver = F
    Strategy.Steps = [solve, correct]
    Step.solve.Operation = SOLVE
    Step.solve.Model.Sources = []
    Step.solve.Model.Beam.Enable = T
    Step.solve.Model.Beam.UseChannelFreq = T
    Step.solve.Model.Cache.Enable = T
    Step.solve.Model.Bandpass.Enable = F
    Step.solve.Model.Phasors.Enable = T
    Step.solve.Model.Gain.Enable = T
    Step.solve.Model.DirectionalGain.Enable = F
    Step.solve.Solve.Mode = COMPLEX
    Step.solve.Solve.Parms = ["Gain:0:0:Phase:*","Gain:1:1:Phase:*"]
    Step.solve.Solve.UVRange = [150, 999999]
    Step.solve.Solve.ExclParms = []
    Step.solve.Solve.CalibrationGroups = []
    Step.solve.Solve.CellSize.Freq = 0
    Step.solve.Solve.CellSize.Time = 4
    Step.solve.Solve.CellChunkSize = 100
    Step.solve.Solve.PropagateSolutions = F
    Step.solve.Solve.Options.MaxIter = 500
    Step.solve.Solve.Options.EpsValue = 1e-8
    Step.solve.Solve.Options.EpsDerivative = 1e-8
    Step.solve.Solve.Options.ColFactor = 1e-9
    Step.solve.Solve.Options.LMFactor = 1.0
    Step.solve.Solve.Options.BalancedEqs = F
    Step.solve.Solve.Options.UseSVD = T
    Step.correct.Operation = CORRECT
    Step.correct.Model.Sources = []
    Step.correct.Model.Phasors.Enable = T
    Step.correct.Model.Gain.Enable = T
    Step.correct.Model.Beam.Enable = T
    Step.correct.Model.Beam.UseChannelFreq = T
    Step.correct.Output.Column = CORRECTED_DATA

The output of this step is a ``CORRECTED_DATA`` column for each of the combined target measurement sets and the instrument tables with the phase solutions. These data products are the final products of this pipeline and once they are checked (with e.g. imaging) the intermediate products are no longer required.

.. _figure-phase:

.. figure:: P23-pipeline_target_SB011SB233_phase_all.png
   :figwidth: 90 %

   Example phase solutions after phase calibration of the target field. Top left: phase solutions for all antennas at 130 MHz with reference to CS001. Top right: phase solutions for the remote station RS305 with reference to CS001. The bottom plots are the same as the top plots but at 165 MHz.



.. _pre-facet-calibration-alterations:

pre-facet calibration alterations
---------------------------------

If preparing data for facet calibration alterations are required to the BBS parsets used in Sections  :ref:`section_calib`, :ref:`sec:parmmap_calibtarget` and :ref:`sec:gsmcalibtarget`. The parsets calibtarget.parset, gsmcal.parset and applyparmdb.parset should be used instead of ``calibtarget_beam.parset``, ``gsmcal_beam.parset`` and ``applyparmdb_beam.parset``. These additional parsets only change in their treatment of the LOFAR beam and are given below.

applyparmdb.parset::

    Strategy.InputColumn = DATA
    Strategy.ChunkSize = 1000
    Strategy.UseSolver = F
    Strategy.Steps = [ correct]

    Step.correct.Operation = CORRECT
    Step.correct.Model.Sources = []
    Step.correct.Model.Cache.Enable  = T
    Step.correct.Model.Clock.Enable = T
    Step.correct.Model.Gain.Enable = T
    Step.correct.Model.CommonRotation.Enable = F
    Step.correct.Model.Beam.Enable = T
    Step.correct.Model.Beam.UseChannelFreq = T
    Step.correct.Output.Column = CORRECTED_DATA

calibtarget.parset::

    Strategy.InputColumn = DATA
    Strategy.ChunkSize   = 200
    Strategy.UseSolver   = F
    Strategy.Steps       = [solve]

    Step.solve.Model.Sources                = [] # all in skymodel
    Step.solve.Model.Cache.Enable           = T
    Step.solve.Model.Phasors.Enable         = F
    Step.solve.Model.DirectionalGain.Enable = F
    Step.solve.Model.Gain.Enable            = T
    Step.solve.Model.Clock.Enable           = T
    Step.solve.Model.TEC.Enable             = F
    Step.solve.Operation                    = SOLVE
    Step.solve.Solve.Parms                  = ["Gain:0:0:*","Gain:1:1:*","Clock:*"]
    Step.solve.Solve.CellSize.Freq          = 0
    Step.solve.Solve.CellSize.Time          = 1
    Step.solve.Solve.CellChunkSize          = 100
    Step.solve.Solve.PropagateSolutions     = T
    Step.solve.Solve.Options.MaxIter        = 1
    Step.solve.Solve.Options.EpsValue       = 1e-9
    Step.solve.Solve.Options.EpsDerivative  = 1e-9
    Step.solve.Solve.Options.ColFactor      = 1e-9
    Step.solve.Solve.Options.LMFactor       = 1.0
    Step.solve.Solve.Options.BalancedEqs    = F
    Step.solve.Solve.Options.UseSVD         = T
    Step.solve.Solve.Mode                   = COMPLEX
    Step.solve.Model.Beam.Enable            = F
    Step.solve.Model.Beam.UseChannelFreq    = F

gsmcal.parset::

    Strategy.InputColumn = DATA
    Strategy.ChunkSize = 50
    Strategy.UseSolver = F
    Strategy.Steps = [solve, correct]
    Step.solve.Operation = SOLVE
    Step.solve.Model.Sources = []
    Step.solve.Model.Beam.Enable = T
    Step.solve.Model.Beam.UseChannelFreq = T
    Step.solve.Model.Beam.Mode = ARRAY_FACTOR
    Step.solve.Model.Cache.Enable = T
    Step.solve.Model.Bandpass.Enable = F
    Step.solve.Model.Phasors.Enable = T
    Step.solve.Model.Gain.Enable = T
    Step.solve.Model.DirectionalGain.Enable = F
    Step.solve.Solve.Mode = COMPLEX
    Step.solve.Solve.Parms = ["Gain:0:0:Phase:*","Gain:1:1:Phase:*"]
    Step.solve.Solve.UVRange = [150, 999999]
    Step.solve.Solve.ExclParms = []
    Step.solve.Solve.CalibrationGroups = []
    Step.solve.Solve.CellSize.Freq = 0
    Step.solve.Solve.CellSize.Time = 4
    Step.solve.Solve.CellChunkSize = 25
    Step.solve.Solve.PropagateSolutions = F
    Step.solve.Solve.Options.MaxIter = 500
    Step.solve.Solve.Options.EpsValue = 1e-8
    Step.solve.Solve.Options.EpsDerivative = 1e-8
    Step.solve.Solve.Options.ColFactor = 1e-9
    Step.solve.Solve.Options.LMFactor = 1.0
    Step.solve.Solve.Options.BalancedEqs = F
    Step.solve.Solve.Options.UseSVD = T
    Step.correct.Operation = CORRECT
    Step.correct.Model.Sources = []
    Step.correct.Model.Phasors.Enable = T
    Step.correct.Model.Gain.Enable = T
    Step.correct.Model.Beam.Enable = F
    Step.correct.Model.Beam.UseChannelFreq = F
    Step.correct.Output.Column = CORRECTED_DATA


