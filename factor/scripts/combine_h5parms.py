#! /usr/bin/env python
"""
Script to combine two h5parms
"""
import argparse
from argparse import RawTextHelpFormatter
from losoto.h5parm import h5parm
import sys
import os


def main(h5parm1, h5parm2, outh5parm, solset1='sol000', solset2='sol000'):
    """
    Combines two h5parms

    Parameters
    ----------
    h5parm1 : str
        Filenames of fast-phase h5parm
    h5parm2 : str
        Filenames of slow-gain h5parm
    outh5parm : str
        Filename of the output h5parm
    solset1 : str, optional
        Name of solset for h5parm1
    solset2 : str, optional
        Name of solset for h5parm2
    """
    h1 = h5parm(h5parm1)
    h2 = h5parm(h5parm2)
    ss1 = h1.getSolset(solset=solset1)
    ss2 = h2.getSolset(solset=solset2)

    # Rename slow-phase soltab before combining to avoid conflict with fast-phase soltab
    soltab = ss2.getSoltab('phase000')
    soltab.rename('phase001')

    if os.path.exists(outh5parm):
        os.remove(outh5parm)
    ho = h5parm(outh5parm, readonly=False)

    sso = ho.makeSolset(solsetName = 'sol000', addTables=False)
    ss1.obj._f_copy_children(sso.obj, recursive=True, overwrite=True)
    ss2.obj._f_copy_children(sso.obj, recursive=True, overwrite=True)

    h1.close()
    h2.close()
    ho.close()


if __name__ == '__main__':
    descriptiontext = "Combine two h5parms.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('h51', help='name of input h5 1')
    parser.add_argument('h52', help='name of input h5 2')
    parser.add_argument('outh5', help='name of the output h5')
    args = parser.parse_args()

    main(args.h51, args.h52, args.outh5)
