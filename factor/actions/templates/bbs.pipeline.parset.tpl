pipeline.steps=[prepare, runbbs, finalize]

prepare.control.kind=recipe
prepare.control.type=executable_args
prepare.control.opts.executable={{ prepare_scriptname }}
prepare.control.opts.mapfiles_in=[{{ input_datamap }}, {{ parmdb_datamap }}]
prepare.control.opts.inputkeys=[inputms, parmdb]
prepare.control.opts.arguments=[inputms, parmdb]
prepare.control.opts.max_per_node={{ n_per_node }}

runbbs.control.kind=recipe
runbbs.control.type=executable_args
runbbs.control.opts.mapfiles_in=[{{ vis_datamap }}, {{ skymodel_datamap }}]
runbbs.control.opts.executable={{ lofarroot }}/bin/calibrate-stand-alone
runbbs.control.opts.inputkeys=[inputms, inputskymodel]
runbbs.control.opts.arguments=[{{ flags }}, inputms, {{ parset }}, inputskymodel]
runbbs.control.opts.max_per_node={{ n_per_node }}

finalize.control.kind=recipe
finalize.control.type=executable_args
finalize.control.opts.executable={{ finalize_scriptname }}
finalize.control.opts.mapfiles_in=[{{ input_datamap }}, {{ parmdb_datamap }}]
finalize.control.opts.inputkeys=[inputms, parmdb]
finalize.control.opts.arguments=[inputms, parmdb]
finalize.control.opts.max_per_node={{ n_per_node }}
