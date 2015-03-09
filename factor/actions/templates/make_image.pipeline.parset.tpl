pipeline.steps=[image]

image.control.kind=recipe
image.control.type=casapy
image.control.opts.mapfiles_in=[{{ vis_datamap }}, {{output_datamap }}]
image.control.opts.inputkeys=[inputms, imagebase]

image.control.opts.arguments=[--nologger, --log2term, --nogui, -c, {{ scriptname }}, inputms, imagebase]
image.control.opts.max_per_node={{ ncpu }}
