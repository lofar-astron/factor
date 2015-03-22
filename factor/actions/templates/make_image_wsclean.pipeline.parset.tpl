pipeline.steps=[wsclean]

wsclean.control.kind=recipe
wsclean.control.type=executable_args
wsclean.control.opts.executable={{ imagerroot }}/bin/wsclean
wsclean.control.opts.mapfiles_in=[{{ vis_datamap }}, {{ output_datamap }}]
wsclean.control.opts.inputkeys=[msin, imagename]
wsclean.control.opts.arguments=[-name, imagename, -size, {{ imsize }}, {{ imsize }}, -niter, {{ niter }}, -threshold, {{ threshold }}, -pol, xx, -weight, briggs, {{ robust }}, -scale, {{ cell }}, -joinchannels, -channelsout, {{ nchannels }}, -mgain, {{ 0.8 }}, msin]
wsclean.control.opts.max_per_node={{ ncpu }}
