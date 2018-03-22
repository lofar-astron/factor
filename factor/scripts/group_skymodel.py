#! /usr/bin/env python
"""
Script to group a sky model
"""
import argparse
from argparse import RawTextHelpFormatter
import lsmtool
import sys
import os


def main(model, skymodel):
    """
    Groups a sky model into a single patch

    Parameters
    ----------
    model : str
        Filename of the input makesourcedb sky model
    skymodel : str
        Filename of the output makesourcedb sky model

    """
    s = lsmtool.load(model)
    s.group('single')

    s.write(skymodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Combine two makesourcedb sky models.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('model', help='name of input makesourcedb sky model')
    parser.add_argument('skymodel', help='name of the output makesourcedb sky model')
    args = parser.parse_args()

    main(args.model, args.skymodel)
