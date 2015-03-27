pipeline.steps=[model]

model.control.kind=recipe
model.control.type=executable_args
model.control.opts.mapfiles_in=[{{ input_datamap }}, {{ output_datamap}}]
model.control.opts.inputkeys=[inputskymodel, outputskymodel]
model.control.opts.executable=/usr/bin/python
model.control.opts.arguments=[{{scriptname}}, inputskymodel, outputskymodel, {{ cal_only }}]
model.control.opts.max_per_node={{ n_per_node }}
