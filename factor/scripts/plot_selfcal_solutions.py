#!/usr/bin/python
"""
Script to plot solutions
"""
from losoto.h5parm import h5parm
from losoto.operations import plot
import argparse
from argparse import RawTextHelpFormatter


def main(h5file, refstat=None):
    """
    Plot solutions vs. time

    Parameters
    ----------
    h5file : str
        Name of solution h5parm file
    refstat : str, optional
        Name of referance station. If None, the first stations is used
    """
    h = h5parm(h5file)
    ss = h.getSolset('sol000')

    sols = ['phase000', 'amplitude001', 'phase001']
    roots = ['tec+scalarphase_', 'slow-amplitude_', 'slow-phase_']
    ncols = [1, 0, 0]
    colors = ['', 'pol', 'pol']
    minmaxes = [[-3.2, 3.2], [0, 0], [-3.2, 3.2]]
    for sol, root, ncol, color, minmax in zip(sols, roots, ncols, colors, minmaxes):
        st = ss.getSoltab(sol)
        if sol == 'phase000':
            # add tec000
            stadd = ['tec000']
        else:
            stadd = ''
        if refstat is not None:
            ref = refstat
        else:
            ref = st.ant[0]
        if sol == 'amplitude001':
            ref = ''
        print('Plotting {} solutions...'.format(root[:-1]))
        plot.run(st, ['time'], axisInTable='ant', axisInCol=color, NColFig=ncol, refAnt=ref,
                 prefix=root, minmax=minmax, soltabsToAdd=stadd)
    h.close()


if __name__ == "__main__":
    descriptiontext = "Plot solutions.\n"
    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('h5file', help="Name of solution h5parm file")
    parser.add_argument('--refstat', help="Name of referance station (default: first)", default=None)

    args = parser.parse_args()
    main(args.h5file, args.refstat)
