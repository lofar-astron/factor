pipeline.steps=[fft]

fft.control.kind=recipe
fft.control.type=executable_args
fft.control.opts.executable={{ imagerroot }}/bin/wsclean
fft.control.opts.mapfiles_in=[{{ vis_datamap }}, {{model_datamap}}]
fft.control.opts.inputkeys=[inputms, imagename]
fft.control.opts.arguments=[-predict, -name, imagename, -smallinversion, -size, {{ imsize }}, {{ imsize }}, -scale, {{ cell_deg }}, -channelsout, {{ nchannels }}, -j, {{ ncpu }}, inputms]
fft.control.opts.max_per_node={{ n_per_node }}
