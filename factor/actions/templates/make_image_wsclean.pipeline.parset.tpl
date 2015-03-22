pipeline.steps=[wsclean]

wsclean.control.kind=recipe
wsclean.control.type=executable_args
wsclean.control.opts.executable={{ imagerroot }}/bin/wsclean
wsclean.control.opts.mapfile_in={{ vis_datamap }}
wsclean.control.opts.mapfile_out={{ output_datamap }}
wsclean.control.opts.inputkey=msin
wsclean.control.opts.outputkey=imagename
wsclean.control.opts.arguments=[-name, imagename, -size, {{ imsize }}, {{ imsize }}, -niter, {{ niter }}, -threshold, {{ threshold }}, -pol, xx, -weight, briggs, {{ robust }}, -scale, {{ cell }}, -joinchannels, -channelsout, {{ nchannels }}, -mgain, {{ 0.8 }}, msin]
wsclean.control.opts.max_per_node={{ ncpu }}
