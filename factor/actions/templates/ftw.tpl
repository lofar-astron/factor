import os
import sys
import numpy
import time

# Set up custom ftw task
os.system('cp {{ task_xml_file }} .')
os.system('cp {{ task_py_file }} .')
os.system('buildmytasks')
execfile('mytasks.py')

ms    = sys.argv[6]
modimage  = sys.argv[7] + '.model'
ntermsi   = numpy.int(sys.argv[8])
wplanes   = numpy.int(sys.argv[9])

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

os.system('touch {{ completed_file }}')

