pipeline.steps=[fft]

fft.control.kind=recipe
fft.control.type=executable_args
fft.control.opts.executable=casa
fft.control.opts.mapfiles_in=[{{ vis_datamap }}, {{model_datamap}}, {{ region_datamap }}]
fft.control.opts.inputkeys=[inputms, modelimg, region]
fft.control.opts.arguments=[--nologger, --log2term, --nogui, -c, {{ scriptname }}, inputms, modelimg, {{ nterms }}, {{ imsize }}]
