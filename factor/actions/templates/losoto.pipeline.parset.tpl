pipeline.steps=[import, losoto, export]

import.control.kind=recipe
import.control.type=executable_args
import.control.opts.executable={{ lofarroot }}/bin/H5parm_importer.py
import.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
import.control.opts.inputkeys=[inputms, h5parm]
import.control.opts.arguments=[-i, instrument, h5parm, inputms]

losoto.control.kind=recipe
losoto.control.type=executable_args
losoto.control.opts.executable={{ lofarroot }}/bin/losoto.py
losoto.control.opts.mapfile_in={{ h5parm_datamap }}
losoto.control.opts.inputkey=h5parm
losoto.control.opts.arguments=[h5parm, {{ parset }}]

export.control.kind=recipe
export.control.type=executable_args
export.control.opts.executable={{ lofarroot }}/bin/H5parm_exporter.py
export.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
export.control.opts.inputkeys=[inputms, h5parm]
export.control.opts.arguments=[--clobber, h5parm, inputms, --instrument=instrument]

calscript.control.kind=recipe
calscript.control.type=executable_args
calscript.control.opts.mapfiles_in=[{{ vis_datamap }}, {{ parmdb_datamap }}, {{ skymodel_datamap }}]
calscript.control.opts.executable={{ lofarroot }}/bin/calibrate-stand-alone
calscript.control.opts.inputkeys=[inputms, inputparmdb, inputskymodel]
calscript.control.opts.arguments=[{{ flags }}, --parmdb, inputparmdb, inputms, {{ parset }}, inputskymodel]
