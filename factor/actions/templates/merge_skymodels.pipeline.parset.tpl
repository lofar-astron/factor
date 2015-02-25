pipeline.steps=[merge]

merge.control.kind=recipe
merge.control.type=executable_args
merge.control.opts.mapfiles_in=[{{ skymodel1_datamap }}, {{ skymodel2_datamap }}]
merge.control.opts.inputkeys=[inputskymodel1, inputskymodel2]
merge.control.opts.mapfile_out={{ output_datamap }}
merge.control.opts.outputkey=outputskymodel
merge.control.opts.executable=/usr/bin/python
merge.control.opts.arguments=[{{mergescriptname}}, inputskymodel1, inputskymodel2, outputskymodel]
