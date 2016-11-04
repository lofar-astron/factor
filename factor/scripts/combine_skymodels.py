#! /usr/bin/env python
"""
Script to combine two makesourcedb sky models
"""
import argparse
from argparse import RawTextHelpFormatter
import lsmtool
import sys
import os


def main(calmodel, facetmodel, skymodel, oldmodel=None):
    """
    Combines makesourcedb sky models

    Parameters
    ----------
    calmodel : str
        Filename of the input calibrator makesourcedb sky model
    facetmodel : str
        Filename of the input full-facet makesourcedb sky model
    skymodel : str
        Filename of the output makesourcedb sky model
    oldmodel : str, optional
        Filename of old model to combine with other two

    """
    s_cal = lsmtool.load(calmodel)
    try:
        s_facet = lsmtool.load(facetmodel)
        s_facet.concatenate(s_cal, matchBy='position', radius='10 arcsec', keep='from2')
    except:
        # If facet skymodel is empty, just set s_facet to s_cal
        s_facet = s_cal

    if oldmodel is not None:
        try:
            s_old = lsmtool.load(oldmodel)
            s_facet.concatenate(s_old, matchBy='position', radius='10 arcsec', keep='from1')
        except:
            pass

    s_facet.write(skymodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Combine two makesourcedb sky models.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('calmodel', help='name of input calibrator makesourcedb sky model')
    parser.add_argument('facetmodel', help='name of input full-facet makesourcedb sky model')
    parser.add_argument('skymodel', help='name of the output makesourcedb sky model')
    args = parser.parse_args()

    main(args.calmodel, args.facetmodel, args.skymodel)
