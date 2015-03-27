pipeline.steps=[calscript]

calscript.control.kind=recipe
calscript.control.type=executable_args
calscript.control.opts.mapfiles_in=[{{ vis_datamap }}, {{ skymodel_datamap }}]
calscript.control.opts.executable={{ lofarroot }}/bin/calibrate-stand-alone
calscript.control.opts.inputkeys=[inputms, inputskymodel]
calscript.control.opts.arguments=[{{ flags }}, inputms, {{ parset }}, inputskymodel]
calscript.control.opts.max_per_node={{ n_per_node }}
