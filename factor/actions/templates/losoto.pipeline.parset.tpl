pipeline.steps=[prepare, import, runlosoto, export, finalize]

prepare.control.kind=recipe
prepare.control.type=executable_args
prepare.control.opts.executable={{ prepare_scriptname }}
prepare.control.opts.mapfiles_in=[{{ input_datamap }}, {{ parmdb_datamap }}]
prepare.control.opts.inputkeys=[inputms, parmdb]
prepare.control.opts.arguments=[inputms, parmdb]
prepare.control.opts.max_per_node={{ n_per_node }}

import.control.kind=recipe
import.control.type=executable_args
import.control.opts.executable={{ h5importer_exec }}
import.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
import.control.opts.inputkeys=[inputms, h5parm]
import.control.opts.arguments=[h5parm, inputms, -i, instrument]
import.control.opts.max_per_node={{ n_per_node }}

runlosoto.control.kind=recipe
runlosoto.control.type=executable_args
runlosoto.control.opts.executable={{ losoto_exec }}
runlosoto.control.opts.mapfile_in={{ h5parm_datamap }}
runlosoto.control.opts.inputkey=h5parm
runlosoto.control.opts.arguments=[h5parm, {{ parset }}]
runlosoto.control.opts.max_per_node={{ n_per_node }}

export.control.kind=recipe
export.control.type=executable_args
export.control.opts.executable={{ h5exporter_exec }}
export.control.opts.mapfiles_in=[{{ input_datamap }}, {{ h5parm_datamap }}]
export.control.opts.inputkeys=[inputms, h5parm]
export.control.opts.arguments=[h5parm, inputms, -i, instrument, -c]
export.control.opts.max_per_node={{ ncpu }}

finalize.control.kind=recipe
finalize.control.type=executable_args
finalize.control.opts.executable={{ finalize_scriptname }}
finalize.control.opts.mapfiles_in=[{{ input_datamap }}, {{ parmdb_datamap }}]
finalize.control.opts.inputkeys=[inputms, parmdb]
finalize.control.opts.arguments=[inputms, parmdb, {{ solset }}]
finalize.control.opts.max_per_node={{ n_per_node }}

