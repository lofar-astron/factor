#!/usr/bin/python
"""
Script to plot solutions
"""
from losoto.h5parm import h5parm
from losoto.operations import plot
import argparse
from argparse import RawTextHelpFormatter


def main(h5file, root=None, refstat=None):
    """
    Plot solutions vs. time

    Parameters
    ----------
    h5file : str
        Name of solution h5parm file
    root : str, optional
        Root name for output plots. If None, the soltype is used
    refstat : str, optional
        Name of referance station. If None, the first stations is used
    """
    h = h5parm(h5file)
    ss = h.getSolset('sol000')

    sols = ['tec000', 'phase000', 'amplitude001', 'phase001']
    ncols = [1, 0, 0, 0]
    colors = ['', '', 'pol', 'pol']
    minmaxes = [[0, 0], [-3.2, 3.2], [0, 0], [-3.2, 3.2]]
    for sol, ncol, color, minmax in zip(sols, ncols, colors, minmaxes):
        st = ss.getSoltab(sol)
        if root is None:
            root = soltype + '_'
        ref = st.ant[0]
        if refstat is not None:
            ref = refstat
        print('Plotting {} solutions...'.format(soltype))
        plot.run(st, ['time'], axisInTable='ant', axisInCol=color, NColFig=ncol, refAnt=ref,
                 prefix=root, minmax=minmax)
    h.close()


if __name__ == "__main__":
    descriptiontext = "Plot solutions.\n"
    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('h5file', help="Name of solution h5parm file")
    parser.add_argument('--root', help="Root name for output plots (default: 'soltype_')", default=None)
    parser.add_argument('--refstat', help="Name of referance station (default: first)", default=None)

    args = parser.parse_args()
    main(args.h5file, args.root, args.refstat)
