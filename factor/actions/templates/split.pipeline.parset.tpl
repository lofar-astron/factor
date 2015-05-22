pipeline.steps=[dpppex]

dpppex.control.kind=recipe
dpppex.control.type=dppp
dpppex.control.opts.mapfile_in={{ input_datamap }}
dpppex.control.opts.inputkey=msin.name
dpppex.control.opts.executable={{ lofarroot }}/bin/NDPPP
dpppex.control.opts.mapfile_out={{ output_datamap }}
dpppex.control.opts.outputkey=msout
dpppex.control.opts.max_per_node={{ n_per_node }}

dpppex.parsetarg.msin.datacolumn={{ columnname }}
dpppex.parsetarg.msin.autoweight = False
dpppex.parsetarg.msout.writefullresflag=False
dpppex.parsetarg.steps=[]
