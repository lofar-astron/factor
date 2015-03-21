pipeline.steps=[awimager]

awimager.control.kind=recipe
awimager.control.type=awimager
awimager.control.opts.mapfile_in={{ vis_datamap }}
awimager.control.opts.mapfile_out={{ output_datamap }}
awimager.control.opts.inputkey=data.ms
awimager.control.opts.outputkey=output.imagename
awimager.control.opts.arguments=[]
awimager.control.opts.max_per_node={{ ncpu }}

awimager.parsetarg.data.uvrange='{{ uvrange }}'
awimager.parsetarg.image.nterms={{ nterms }}
awimager.parsetarg.clean.niter={{ niter }}
awimager.parsetarg.clean.gain=0.01
awimager.parsetarg.clean.threshold='{{ threshold }}'
awimager.parsetarg.clean.interactive=False
awimager.parsetarg.clean.imsize=[{{ imsize }}, {{ imsize }}]
awimager.parsetarg.clean.cell=['{{ cell }}', '{{ cell }}']
awimager.parsetarg.weight.type='briggs'
awimager.parsetarg.weight.robust=-0.25
awimager.parsetarg.clean.uvtaper=False
awimager.parsetarg.clean.pbcor=False
awimager.parsetarg.clean.minpb=0.2
awimager.parsetarg.clean.nscales={{ nscales }}
awimager.parsetarg.clean.uservector={{ scales }}
awimager.parsetarg.clean.mask={{ mask }}
