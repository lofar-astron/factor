import sys
try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess
import os

model = sys.argv[1]
nterms = sys.argv[2]
outfile = sys.argv[3]

cmd = 'casapy2bbs2.py -t {0} {1} {2}'.format(nterms, model, outfile)
p = subprocess.Popen(cmd, shell=True)
p.wait()
