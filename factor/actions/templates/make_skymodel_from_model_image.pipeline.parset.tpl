pipeline.steps=[model]

model.control.kind=recipe
model.control.type=executable_args
model.control.opts.executable=/usr/bin/python
model.control.opts.mapfile_in={{ input_datamap }}
model.control.opts.inputkey=inputmodel
model.control.opts.arguments=[{{ scriptname }}, inputmodel, {{ nterms }}]
