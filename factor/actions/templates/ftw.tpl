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
imsizep   = numpy.int(sys.argv[9])

wplanes   = 1

if imsizep > 512:
    wplanes = 64

if imsizep > 799:
    wplanes = 96

if imsizep > 1023:
    wplanes = 128

if imsizep > 1599:
    wplanes = 256

if imsizep > 2047:
    wplanes = 384

if imsizep > 3000:
    wplanes = 448

if imsizep > 4095:
    wplanes = 512

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

