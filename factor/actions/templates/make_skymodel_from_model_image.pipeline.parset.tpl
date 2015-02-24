pipeline.steps=[model]

model.control.kind=recipe
model.control.type=executable_args
model.control.opts.executable=casapy2bbs2.py
model.control.opts.mapfiles_in=[{{ input_datamap_model_images }}, {{ input_datamap_mask_images}}]
model.control.opts.inputkeys=[inputmodel, inputmask]
model.control.opts.mapfile_out={{ output_datamap }}
model.control.opts.outputkey=outputmodel
model.control.opts.arguments=[-t {{ nterms }}, -m {{ inputmask }}, inputmodel, outputmodel]
