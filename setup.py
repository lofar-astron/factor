from __future__ import print_function
from setuptools import setup, Command
import os
import factor._version

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
    version=factor._version.__version__,
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
    install_requires=['numpy', 'scipy', 'astropy', 'jinja2', 'aplpy>=1.0', 'LSMTool>=1.1', ],
    dependency_links=['https://github.com/darafferty/LSMTool'],
    scripts = ['bin/runfactor','bin/checkfactor'],
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
