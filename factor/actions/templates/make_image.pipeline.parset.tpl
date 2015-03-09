pipeline.steps=[image]

image.control.kind=recipe
image.control.type=casapy
image.control.opts.mapfile_in={{ vis_datamap }}
image.control.opts.inputkey=inputms
image.control.opts.arguments=[--nologger, --log2term, --nogui, -c, {{ scriptname }}, inputms]
image.control.opts.max_per_node={{ ncpu }}
