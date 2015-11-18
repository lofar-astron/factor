"""
Script to perform FT (to be run only by casapy)

The assumed call is:

    casapy --nologger --nogui -c do_ft.py  modelimg nterms wplanes task_xml_file task_py_file inputmss

modelimg      : CASA model image. If nterms > 1, this should be the root name
                (without the .tt0, .tt1, etc.)
nterms        : nterms in clean
wplanes       : number of wplanes
task_xml_file : full path to the ftw.xml file
task_py_file  : full path to the task_ftw.py file
inputmss      : MS files receiving model vis

"""
import os
import sys
import numpy
import time

if len(sys.argv) < 11:
    print "do_ft_multi.py: Too few arguments!"

# This is currently a constant
files_to_concat = 100

# Set up custom ftw task
os.system('cp {0} .'.format(sys.argv[8]))
os.system('cp {0} .'.format(sys.argv[9]))
if 'INSTALLDIR' in os.environ:
    # Remove this env variable, otherwise CASA puts "mytasks.py" there
    del os.environ['INSTALLDIR']
os.system('buildmytasks')
execfile('mytasks.py')

modimage  = sys.argv[5] + '.model'
ntermsi   = numpy.int(sys.argv[6])
wplanes   = numpy.int(sys.argv[7])

mslist = []
for argstring in sys.argv[10:]:
    mslist.extend([ms.strip(" []\'\"") for ms in argstring.split(',')])

if ntermsi == 1:
    mod = [modimage]
if ntermsi == 2:
    mod = [modimage+'.tt0', modimage+'.tt1']
if ntermsi == 3:
    mod = [modimage+'.tt0', modimage+'.tt1', modimage+'.tt2']

for ms in mslist:    
    if wplanes > 1:
        ftw(vis=ms, field="", spw="", model=mod, nterms=ntermsi, reffreq="",
            wprojplanes=wplanes, complist="", incremental=False, usescratch=True,
            async=False)
    else:
        ft(vis=ms, field="", spw="", model=mod, nterms=ntermsi, reffreq="",
           complist="", incremental=False, usescratch=True, async=False)
