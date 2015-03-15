pipeline.steps=[mask]

mask.control.kind=recipe
mask.control.type=executable_args
mask.control.opts.executable=/usr/bin/python
mask.control.opts.mapfile_in={{ input_datamap }}
mask.control.opts.inputkey=inputms
mask.control.opts.arguments=[{{ maskscriptname }}, inputms]
