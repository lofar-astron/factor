import os
import sys
import numpy
import time

# Set up custom ftw task
!cp {{ task_xml_file }} .
!cp {{ task_py_file }} .
!buildmytasks
execfile('mytasks.py')

ms    = sys.argv[4]
modimage  = sys.argv[5]
ntermsi   = numpy.int(6)
imsizep   = numpy.int(7)
if len(sys.argv) > 8:
    region = sys.argv[8]
else:
    region = None

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

if region is not None:
    # Trim model image to region
    pass

if wplanes > 1:
    ftw(vis=ms, field="", spw="", model=mod, nterms=ntermsi, reffreq="",
    wprojplanes=wplanes, complist="", incremental=False, usescratch=True,
    async=False)
else:
    ft(vis=ms, field="", spw="", model=mod, nterms=ntermsi, reffreq="",
    complist="", incremental=False, usescratch=True, async=False)

