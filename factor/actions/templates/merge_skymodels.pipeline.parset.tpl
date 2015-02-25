pipeline.steps=[merge]

merge.control.kind=recipe
merge.control.type=executable_args
merge.control.opts.executable=/usr/bin/python
merge.control.opts.mapfiles_in=[{{ skymodel1_datamap }}, {{ skymodel2_datamap }}, {{ output_datamap }}]
merge.control.opts.inputkeys=[inputskymodel1, inputskymodel2, outputskymodel]
merge.control.opts.arguments=[{{mergescriptname}}, inputskymodel1, inputskymodel2, outputskymodel]
