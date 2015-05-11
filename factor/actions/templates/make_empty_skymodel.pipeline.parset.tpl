pipeline.steps=[model]

model.control.kind=recipe
model.control.type=executable_args
model.control.opts.mapfiles_in=[{{ model_datamap }}]
model.control.opts.inputkeys=[skymodel]
model.control.opts.executable=/usr/bin/touch
model.control.opts.arguments=[skymodel]
model.control.opts.max_per_node={{ n_per_node }}
