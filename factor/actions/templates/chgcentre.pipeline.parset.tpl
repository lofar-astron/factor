pipeline.steps=[chgcentre]

chgcentre.control.kind=recipe
chgcentre.control.type=executable_args
chgcentre.control.opts.executable={{ imagerroot }}/bin/chgcentre
chgcentre.control.opts.mapfile_in={{ output_datamap }}
chgcentre.control.opts.inputkey=inputms
chgcentre.control.opts.arguments=[{{ scriptname }}, -zenith, -shiftback, inputms]
