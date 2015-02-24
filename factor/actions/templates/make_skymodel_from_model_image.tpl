import sys
import subprocess

model = sys.argv[1]
nterms = sys.argv[2]
skymodel = model.split('.model')[0] + '.skymodel'

cmd = 'casapy2bbs2.py -t {0} {1} {2}'.format(nterms, model, skymodel)
p = subprocess.Popen(cmd, shell=True)
p.wait()
