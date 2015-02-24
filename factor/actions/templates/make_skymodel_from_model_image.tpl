import sys
import subprocess
import os

model = sys.argv[1]
nterms = sys.argv[2]
outputdir = sys.argv[3]
skymodel = os.path.basename(model.split('.model')[0]) + '.skymodel'
outfile = os.path.join(outputdir, skymodel)

cmd = 'casapy2bbs2.py -t {0} {1} {2}'.format(nterms, model, outfile)
p = subprocess.Popen(cmd, shell=True)
p.wait()
