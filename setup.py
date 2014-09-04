from __future__ import print_function
from setuptools import setup, Command
import os
import factor._version


description = 'Facet calibration for LOFAR'
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
    name='factor',
    version=factor._version.__version__,
    url='http://github.com/revoltek/facotr/',
    author='Francesco de Gasperin',
    author_email='fdg@hs.uni-hamburg.de',
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
    install_requires=['numpy'],
    scripts = ['bin/factor.py'],
    packages=['facotor','factor.operations'],
    cmdclass = {'test': PyTest},
    )
