pipeline.steps=[fft]

fft.control.kind=recipe
fft.control.type=casapy
fft.control.opts.mapfiles_in=[{{ vis_datamap }}, {{model_datamap}}, {{ completed_datamap }}]
fft.control.opts.inputkeys=[inputms, modelimg]
fft.control.opts.arguments=[--nologger, --log2term, --nogui, -c, {{ scriptname }}, inputms, modelimg, completedfile, {{ nterms }}, {{ wplanes }}]
fft.control.opts.max_per_node={{ n_per_node }}
