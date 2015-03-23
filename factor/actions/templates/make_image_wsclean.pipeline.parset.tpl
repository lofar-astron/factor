pipeline.steps=[wsclean]

wsclean.control.kind=recipe
wsclean.control.type=executable_args
wsclean.control.opts.executable={{ imagerroot }}/bin/wsclean
wsclean.control.opts.mapfiles_in=[{{ vis_datamap }}, {{ output_datamap }}]
wsclean.control.opts.inputkeys=[msin, imagename]
wsclean.control.opts.arguments=[-name, imagename, -size, {{ imsize }}, {{ imsize }}, -niter, {{ niter }}, -threshold, {{ threshold_jy }}, -weight, briggs, -0.25, -scale, {{ cell_deg }}, -joinchannels, -channelsout, {{ nchannels }}, -mgain, {{ 0.8 }}, {{wsclean_multiscale}} {{mask}} msin]
wsclean.control.opts.max_per_node={{ ncpu }}
