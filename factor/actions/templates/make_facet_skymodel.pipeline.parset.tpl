pipeline.steps=[model]

model.control.kind=recipe
model.control.type=executable_args
model.control.opts.mapfile_in={{ input_datamap }}
model.control.opts.inputkey=inputskymodel
model.control.opts.mapfile_out={{ output_datamap }}
model.control.opts.outputkey=outputskymodel
model.control.opts.executable=/usr/bin/python
model.control.opts.arguments=[{{scriptname}}, inputskymodel, outputskymodel, {{ cal_only }}]

