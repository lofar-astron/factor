pipeline.steps=[h5import, losoto, h5export]

h5import.control.kind=recipe
h5import.control.type=executable_args
h5import.control.opts.executable={{ lofarroot }}/bin/H5parm_importer.py
h5import.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
h5import.control.opts.inputkeys=[inputms, h5parm]
h5import.control.opts.arguments=[--instrument=instrument, h5parm, inputms]

losoto.control.kind=recipe
losoto.control.type=executable_args
losoto.control.opts.executable={{ lofarroot }}/bin/losoto.py
losoto.control.opts.mapfile_in={{ h5parm_datamap }}
losoto.control.opts.inputkey=h5parm
losoto.control.opts.arguments=[h5parm, {{ parset }}]

h5export.control.kind=recipe
h5export.control.type=executable_args
h5export.control.opts.executable={{ lofarroot }}/bin/H5parm_exporter.py
h5export.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
h5export.control.opts.inputkeys=[inputms, h5parm]
h5export.control.opts.arguments=[--clobber, h5parm, inputms, --instrument=instrument]
