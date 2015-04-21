pipeline.steps=[casapy]

casapy.control.kind=recipe
casapy.control.type=casapy
casapy.control.opts.mapfile_in={{ vis_datamap }}
casapy.control.opts.mapfile_out={{ output_datamap }}
casapy.control.opts.inputkey=inputms
casapy.control.opts.outputkey=outputimage
casapy.control.opts.arguments=[--nologger, --log2term, --nogui, -c, {{ scriptname }}, inputms, outputimage]
casapy.control.opts.max_per_node={{ n_per_node }}
