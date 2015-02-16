[DEFAULT]
lofarroot = {{ lofarroot }}
pythonpath = {{ lofarpythonpath }}
runtime_directory = {{ runtime_dir }}
recipe_directories = [%(pythonpath)s/lofarpipe/recipes, /home/sttf201/PipelineExample]
working_directory = {{ working_dir }}
task_files = [%(lofarroot)s/share/pipeline/tasks.cfg, /home/sttf201/PipelineExample/tasks.cfg]

[layout]
job_directory = %(runtime_directory)s/%(job_name)s

[cluster]
clusterdesc = %(lofarroot)s/share/local.clusterdesc

[deploy]
engine_ppath = %(pythonpath)s:%(pyraproot)s/lib:/opt/cep/pythonlibs/lib/python/site-packages
engine_lpath = %(lofarroot)s/lib:%(casaroot)s/lib:%(pyraproot)s/lib:%(hdf5root)s/lib:%(wcsroot)s/lib

[logging]
log_file = %(runtime_directory)s/%(job_name)s/logs/%(start_time)s/pipeline.log
xml_stat_file = %(runtime_directory)s/%(job_name)s/logs/%(start_time)s/statistics.xml

[remote]
method = local
max_per_node = {{ ncpu }}

