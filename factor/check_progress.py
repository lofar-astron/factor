"""
Module that checks the progress of a run
"""
import sys
import os
import shutil
import numpy as np
import logging
import factor
import factor.directions
import factor.parset
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon


log = logging.getLogger('progress')


def run(parset_file):
    """
    Displays the current progress of a run
    """
    directions_list = load_directions(parset_file)

    plot_state(directions_list)


def load_directions(parset_file):
    """
    Return directions for a run
    """
    # Read parset
    parset = factor.parset.parset_read(parset_file)

    # Load directions. First check for user-supplied directions file then for
    # Factor-generated file from a previous run
    dir_parset = parset['direction_specific']
    if 'directions_file' in dir_parset:
        directions = factor.directions.directions_read(dir_parset['directions_file'],
            parset['dir_working'])
    elif os.path.exists(os.path.join(parset['dir_working'], 'factor_directions.txt')):
        directions = factor.directions.directions_read(os.path.join(parset['dir_working'],
            'factor_directions.txt'), parset['dir_working'])
    for direction in directions:
        direction.load_state()

    return directions


def plot_state(directions_list):
    """
    Plots the facets of a run
    """
    # Set up coordinate system and figure
    points, midRA, midDec = factor.directions.getxy(directions_list)
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot facets
    for direction in directions_list:
        RAverts = direction.vertices[0]
        Decverts = direction.vertices[1]
        xverts, yverts = factor.directions.radec2xy(RAverts, Decverts,
            refRA=midRA, refDec=midDec)
        xyverts = [np.array([xp, yp]) for xp, yp in zip(xverts, yverts)]
        mpl_poly = Polygon(np.array(poly), facecolor="g", lw=0, alpha=0.4)
        ax.add_patch(mpl_poly)

    ax.relim()
    ax.autoscale()
    plt.show()


#         selfcal_images = find_selfcal_image(direction)
