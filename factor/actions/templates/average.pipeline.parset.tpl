pipeline.steps=[dpppex]

dpppex.control.kind=recipe
dpppex.control.type=dppp
dpppex.control.opts.mapfile_in={{ input_datamap }}
dpppex.control.opts.inputkey=msin
dpppex.control.opts.executable={{ lofarroot }}/bin/NDPPP
dpppex.control.opts.mapfile_out={{ output_datamap }}
dpppex.control.opts.outputkey=msout
dpppex.control.opts.max_per_node={{ n_per_node }}

dpppex.parsetarg.msin.datacolumn={{ columnname }}
dpppex.parsetarg.steps=[flag, avg]
dpppex.parsetarg.flag.type=aoflagger
dpppex.parsetarg.avg.type=squash
dpppex.parsetarg.avg.freqstep={{ freqstep }}
dpppex.parsetarg.avg.timestep={{ timestep }}
