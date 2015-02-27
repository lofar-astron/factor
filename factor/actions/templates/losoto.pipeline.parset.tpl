pipeline.steps=[import, losoto, export]

import.control.kind=recipe
import.control.type=executable_args
import.control.opts.executable=/usr/bin/python
import.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
import.control.opts.inputkeys=[inputms, h5parm]
import.control.opts.arguments=[{{ lofarroot }}/bin/H5parm_importer.py, h5parm, inputms, --instrument=instrument]

losoto.control.kind=recipe
losoto.control.type=executable_args
losoto.control.opts.executable=/usr/bin/python
losoto.control.opts.mapfile_in={{ h5parm_datamap }}
losoto.control.opts.inputkey=h5parm
losoto.control.opts.arguments=[{{ lofarroot }}/bin/losoto.py, h5parm, {{ parset }}]

export.control.kind=recipe
export.control.type=executable_args
export.control.opts.executable=/usr/bin/python
export.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
export.control.opts.inputkeys=[inputms, h5parm]
export.control.opts.arguments=[{{ lofarroot }}/bin/H5parm_exporter.py, --clobber, h5parm, inputms, --instrument=instrument]
