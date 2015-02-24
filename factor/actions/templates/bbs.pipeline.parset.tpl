pipeline.steps=[calscript]

calscript.control.kind=recipe
calscript.control.type=executable_args
calscript.control.opts.mapfiles_in=[{{ vis_datamap }}, {{ parmdb_datamap }}, {{ skymodel_datamap }}]
calscript.control.opts.executable={{ lofarroot }}/bin/calibrate-stand-alone
calscript.control.opts.inputkeys=[inputms, inputparmdb, inputskymodel]
calscript.control.opts.arguments=[{{ flags }}, --parmdb=inputparmdb, inputms, {{ parset }}, inputskymodel]
