"""
Script to perform FT (to be run only by casapy)

The assumed call is:

    casapy --nologger -c ftw.py inputms modelimg nterms wplanes task_xml_file task_py_file

inputms : MS file receiving model vis
modelimg : CASA model image. If nterms > 1, this should be the root name
    (without the .tt0, .tt1, etc.)
nterms : nterms in clean
wplanes : number of wplanes
task_xml_file : full path to the ftw.xml file
task_py_file : full path to the task_ftw.py file

Note that the script is wrapped in a try/except block to ensure that casapy
returns a proper return code and that we need to use os._exit(), not sys.exit()
or exit().

"""
import os
import sys
import numpy
import time

try:
    # Set up custom ftw task
    os.system('cp {0} .'.format(sys.argv[8]))
    os.system('cp {0} .'.format(sys.argv[9]))
    os.system('buildmytasks')
    execfile('mytasks.py')

    ms    = sys.argv[4]
    modimage  = sys.argv[5] + '.model'
    ntermsi   = numpy.int(sys.argv[6])
    wplanes   = numpy.int(sys.argv[7])

    if ntermsi == 1:
        mod = [modimage]
    if ntermsi == 2:
        mod = [modimage+'.tt0', modimage+'.tt1']
    if ntermsi == 3:
        mod = [modimage+'.tt0', modimage+'.tt1', modimage+'.tt2']

    if wplanes > 1:
        ftw(vis=ms, field="", spw="", model=mod, nterms=ntermsi, reffreq="",
        wprojplanes=wplanes, complist="", incremental=False, usescratch=True,
        async=False)
    else:
        ft(vis=ms, field="", spw="", model=mod, nterms=ntermsi, reffreq="",
        complist="", incremental=False, usescratch=True, async=False)
except SystemExit:
    pass
except:
    os._exit(1)
