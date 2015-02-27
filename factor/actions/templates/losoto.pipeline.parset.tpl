pipeline.steps=[import, losoto, export]

import.control.kind=recipe
import.control.type=executable_args
import.control.opts.executable=/usr/bin/python
import.control.opts.mapfile_in={{ input_datamap }}
import.control.opts.inputkey=inputms
import.control.opts.mapfile_out={{ h5parm_datamap }}
import.control.opts.inputkey=h5parm
import.control.opts.arguments=[{{ lofarroot }}/bin/H5parm_importer.py, h5parm, inputms]

losoto.control.kind=recipe
losoto.control.type=executable_args
losoto.control.opts.executable=/usr/bin/python
losoto.control.opts.mapfile_in={{ input_datamap }}
losoto.control.opts.inputkey=h5parm
losoto.control.opts.arguments=[{{ lofarroot }}/bin/losoto.py, h5parm, {{ parset }}]

export.control.kind=recipe
export.control.type=executable_args
export.control.opts.executable=/usr/bin/python
export.control.opts.mapfile_in={{ input_datamap }}
export.control.opts.inputkey=inputms
export.control.opts.mapfile_out={{ h5parm_datamap }}
export.control.opts.outputkey=h5parm
export.control.opts.arguments=[{{ lofarroot }}/bin/H5parm_exporter.py, -c, h5parm, inputms, -i, instrument]

