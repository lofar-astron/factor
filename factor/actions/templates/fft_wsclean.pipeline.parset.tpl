pipeline.steps=[fft]

fft.control.kind=recipe
fft.control.type=executable_args
fft.control.opts.executable={{ imagerroot }}/bin/wsclean
fft.control.opts.mapfiles_in=[{{ vis_datamap }}, {{model_datamap}}]
fft.control.opts.inputkeys=[inputms, imagename]
fft.control.opts.arguments=[-predict, -name, imagename, -size, {{ imsize }}, {{ imsize }}, -niter, {{ niter }}, -threshold, {{ threshold }}, -pol, xx, -weight, briggs, {{ robust }}, -scale, {{ cell }}, -joinchannels, -channelsout, {{ nchannels }}, -mgain, {{ 0.8 }}, inputms]
fft.control.opts.max_per_node={{ ncpu }}
