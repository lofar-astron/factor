pipeline.steps=[pad, fft]

pad.control.kind=recipe
pad.control.type=executable_args
pad.control.opts.executable={{ padscriptname }}
pad.control.opts.mapfiles_in=[{{ model_datamap }}]
pad.control.opts.inputkeys=[inputimage]
pad.control.opts.arguments=[inputimage]
pad.control.opts.max_per_node={{ n_per_node }}

fft.control.kind=recipe
fft.control.type=executable_args
fft.control.opts.executable={{ imagerroot }}/bin/wsclean
fft.control.opts.mapfiles_in=[{{ vis_datamap }}, {{ padded_model_datamap }}]
fft.control.opts.inputkeys=[inputms, imagename]
fft.control.opts.arguments=[-predict, -name, imagename, -size, {{ imsize }}, {{ imsize }}, -scale, {{ cell_deg }}, -channelsout, {{ nchannels }}, -j, {{ ncpu }}, inputms]
fft.control.opts.max_per_node={{ n_per_node }}
