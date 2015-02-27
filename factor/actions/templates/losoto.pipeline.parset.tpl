pipeline.steps=[import, losoto, export]

import.control.kind=recipe
import.control.type=executable_args
import.control.opts.executable={{ lofarroot }}/bin/H5parm_importer.py
import.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
import.control.opts.inputkeys=[inputms, h5parm]
import.control.opts.arguments=[h5parm, inputms, --instrument=instrument]
import.control.opts.parsetasfile=True

losoto.control.kind=recipe
losoto.control.type=executable_args
losoto.control.opts.executable={{ lofarroot }}/bin/losoto.py
losoto.control.opts.mapfile_in={{ h5parm_datamap }}
losoto.control.opts.inputkey=h5parm
losoto.control.opts.arguments=[h5parm, {{ parset }}]
losoto.control.opts.parsetasfile=True

export.control.kind=recipe
export.control.type=executable_args
export.control.opts.executable={{ lofarroot }}/bin/H5parm_exporter.py
export.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
export.control.opts.inputkeys=[inputms, h5parm]
export.control.opts.arguments=[--clobber, h5parm, inputms, --instrument=instrument]
export.control.opts.parsetasfile=True
