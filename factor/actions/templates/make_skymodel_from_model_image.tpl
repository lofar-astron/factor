import sys

model = sys.argv[1]
nterms = sys.argv[2]
skymodel = model.split('.model')[0] + '.skymodel'

casapy2bbs2.py -t {0} {1} {2}'.format(nterms, model, skymodel)
