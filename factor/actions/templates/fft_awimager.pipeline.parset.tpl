pipeline.steps=[fft]

fft.control.kind=recipe
fft.control.type=awimager
fft.control.opts.mapfiles_in=[{{ vis_datamap }}, {{model_datamap}}]
fft.control.opts.inputkeys=[data.ms, predict.modelimage]
fft.control.opts.arguments=[]
fft.control.opts.max_per_node={{ ncpu }}

awimager.parsetarg.operation = predict
