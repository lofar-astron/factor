#! /usr/bin/env python
"""
Script to combine two makesourcedb sky models
"""
import argparse
from argparse import RawTextHelpFormatter
import lsmtool
import sys
import os


def main(model1, model2, skymodel):
    """
    Combines makesourcedb sky models

    Parameters
    ----------
    model1 : str
        Filename of the input makesourcedb sky model 1
    model2 : str
        Filename of the input makesourcedb sky model 2
    skymodel : str
        Filename of the output makesourcedb sky model

    """
    try:
        s1 = lsmtool.load(model1)
    except:
        # If first sky model is empty or cannot be loaded, just copy second one
        # to output file
        os.system('cp -f {0} {1}'.format(model2, skymodel))
        return

    # Now try to load second sky model and combine with first one
    try:
        s2 = lsmtool.load(model2)

        # Combine sky models, keeping all sources
        s1.concatenate(s2, keep='all')
    except:
        # If second sky model is empty or cannot be loaded, just save s1 to output
        pass

    s1.write(skymodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Combine two makesourcedb sky models.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('model1', help='name of input makesourcedb sky model 1')
    parser.add_argument('model2', help='name of input makesourcedb sky model 2')
    parser.add_argument('skymodel', help='name of the output makesourcedb sky model')
    args = parser.parse_args()

    main(args.model1, args.model2, args.skymodel)
