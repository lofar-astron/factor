pipeline.steps=[awimager]

awimager.control.kind=recipe
awimager.control.type=awimager
awimager.control.opts.mapfile_in={{ vis_datamap }}
awimager.control.opts.mapfile_out={{ output_datamap }}
awimager.control.opts.inputkey=data.ms
awimager.control.opts.outputkey=output.imagename
awimager.control.opts.arguments=[]
awimager.control.opts.max_per_node={{ ncpu }}

awimager.parsetarg.operation = clean
awimager.parsetarg.data.uvrange='{{ uvrange }}'
awimager.parsetarg.image.nterms={{ nterms }}
awimager.parsetarg.image.cellsize='{{ cell }}'
awimager.parsetarg.image.npix={{ imsize }}
awimager.parsetarg.weight.type='robust'
awimager.parsetarg.weight.robust=-0.25
awimager.parsetarg.clean.niter={{ niter }}
awimager.parsetarg.clean.gain=0.01
awimager.parsetarg.clean.threshold='{{ threshold }}'
awimager.parsetarg.clean.nscales={{ nscales }}
awimager.parsetarg.clean.uservector={{ scales }}
awimager.parsetarg.clean.maskimage='{{ mask }}'
