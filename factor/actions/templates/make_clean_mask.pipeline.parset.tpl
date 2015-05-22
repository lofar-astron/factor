pipeline.steps=[mask]

mask.control.kind=recipe
mask.control.type=executable_args
mask.control.opts.executable=/usr/bin/python
mask.control.opts.mapfiles_in=[{{ input_datamap }}, {{ output_datamap }}]
mask.control.opts.inputkeys=[inputms, outputmask]
mask.control.opts.arguments=[{{ scriptname }}, inputms, outputmask]
