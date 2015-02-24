mask.control.kind=recipe
mask.control.type=executable_args
mask.control.opts.executable=casapy2bbs2.py
mask.control.opts.mapfiles_in=[{{ input_datamap_model_images }}, {{ input_datamap_mask_images}}]
mask.control.opts.inputkeys=[inputmodel, inputmask]
mask.control.opts.mapfile_out={{ output_datamap }}
mask.control.opts.outputkey=outputmodel
mask.control.opts.arguments=[-t {{ nterms }}, -m {{ inputmask }}, inputmodel, outputmodel]
