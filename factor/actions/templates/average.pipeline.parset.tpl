pipeline.steps=[dpppex]

dpppex.control.kind=recipe
dpppex.control.type=dppp
dpppex.control.opts.mapfile_in=createmap.output.mapfile
dpppex.control.opts.executable={{ lofarroot }}/bin/NDPPP
dpppex.control.opts.inputkey=msin
dpppex.control.opts.outputkey=msout

dpppex.parsetarg.msin.datacolumn={{ columnname }}
dpppex.parsetarg.steps=[avg]
dpppex.parsetarg.avg.type=squash
dpppex.parsetarg.avg.freqstep={{ freqstep }}
dpppex.parsetarg.avg.timestep={{ timestep }}

