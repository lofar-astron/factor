import os
import glob

# automatically import all the other file in this directory
__all__ = [ os.path.basename(f)[:-3] for f in glob.glob(os.path.dirname(__file__)+"/*.py") if not f.endswith('__init__.py')]
for x in __all__: __import__(x, locals(), globals())
