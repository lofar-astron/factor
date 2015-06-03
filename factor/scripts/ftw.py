"""
Script to perform FT (to be run only by casapy)

The assumed call is:

casapy --nologger -c ftw.py inputms modelimg nterms wplanes task_xml_file task_py_file

where:

Parameters
----------
inputms : str
    MS file receiving model vis
modelimg : str
    CASA model image. If nterms > 1, this should be the root name (without the
    .tt0, .tt1, etc.)
nterms : int
wplanes : int
task_xml_file : str
    Full path to the ftw.xml file
task_py_file : str
    Full path to the task_ftw.py file

"""
import os
import sys
import numpy
import time

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
