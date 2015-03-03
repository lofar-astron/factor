pipeline.steps=[calscript]

calscript.control.kind=recipe
calscript.control.type=executable_args
calscript.control.opts.mapfile_in={{ vis_datamap }}
calscript.control.opts.executable={{ lofarroot }}/bin/calibrate-stand-alone
calscript.control.opts.inputkey=inputms
calscript.control.opts.arguments=[{{ flags }}, inputms, {{ parset }}]]
