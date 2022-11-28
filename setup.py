from __future__ import print_function
from setuptools import setup, Command
import os

# Functions read() and get_version() were copied from Pip package.
# Purpose is to get version info from current package without it
# being installed (which is usually the case when setup.py is run).
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            # __version__ = "0.9"
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")

description = 'FACTOR: Facet calibration for LOFAR'
long_description = description
if os.path.exists('README.md'):
    with open('README.md') as f:
        long_description=f.read()


class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import sys,subprocess
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)


setup(
    name='FACTOR',
    version=get_version('factor/_version.py'),
    url='http://github.com/lofar-astron/factor/',
    description=description,
    long_description=long_description,
    platforms='any',
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: 1 - Alpha',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
    install_requires=['numpy', 'scipy', 'astropy', 'jinja2', 'aplpy>=1.0', 'LSMTool>=1.2'],
    dependency_links=['https://github.com/darafferty/LSMTool'],
    scripts = ['bin/runfactor','bin/checkfactor','bin/archivefactor','bin/unarchivefactor'],
    packages=['factor', 'factor.operations', 'factor.lib'],
    package_data={'factor': [
        'parsets/*',
        'pipeline/*.cfg',
        'pipeline/parsets/*',
        'pipeline/plugins/*',
        'pipeline/recipes/nodes/*',
        'scripts/*',
        'skymodels/*']},
    cmdclass = {'test': PyTest},
    )
