pipeline.steps=[casapy]

casapy.control.kind=recipe
casapy.control.type=casapy
casapy.control.opts.mapfiles_in=[{{ vis_datamap }}, {{ output_datamap }}]
casapy.control.opts.inputkeys=[inputms, outputimage]
casapy.control.opts.arguments=[--nologger, --log2term, --nogui, -c, {{ scriptname }}, inputms, outputimage]
casapy.control.opts.max_per_node={{ n_per_node }}

