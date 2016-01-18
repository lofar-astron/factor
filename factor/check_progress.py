"""
Module that checks the progress of a run
"""
import sys
import os
import shutil
import numpy as np
import pickle
import logging
import glob
import ConfigParser
import factor
import factor.directions
import factor.parset
import pyrap.images as pim
import ast
from lofarpipe.support.data_map import DataMap
try:
    from matplotlib import pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.ticker import FuncFormatter
except Exception as e:
    raise ImportError('PyPlot could not be imported. Plotting is not '
        'available: {0}'.format(e.message))
try:
    from wcsaxes import WCSAxes
    hasWCSaxes = True
except:
    hasWCSaxes = False

log = logging.getLogger('progress')


def run(parset_file):
    """
    Displays the current progress of a run
    """
    directions_list = load_directions(parset_file)

    log.info('Plotting facets...')
    log.info('Left-click on a facet to see its current state')
    log.info('Right-click on a facet to display its images\n')
    plot_state(directions_list)


def load_directions(parset_file):
    """
    Return directions for a run
    """
    # Read parset
    parset = parset_read(parset_file)

    # Load directions. First check for user-supplied directions file then for
    # Factor-generated file from a previous run
    direction_list = []
    dir_parset = parset['direction_specific']
    if 'directions_file' in dir_parset:
        directions = factor.directions.directions_read(dir_parset['directions_file'],
            parset['dir_working'])
    elif os.path.exists(os.path.join(parset['dir_working'], 'factor_directions.txt')):
        directions = factor.directions.directions_read(os.path.join(parset['dir_working'],
            'factor_directions.txt'), parset['dir_working'])
    for direction in directions:
        has_state = direction.load_state()
        if has_state:
            direction_list.append(direction)

    return direction_list


def plot_state(directions_list):
    """
    Plots the facets of a run
    """
    global midRA, midDec

    # Set up coordinate system and figure
    points, midRA, midDec = factor.directions.getxy(directions_list)
    fig = plt.figure(1, figsize=(7.66,7))
    if hasWCSaxes:
        wcs = factor.directions.makeWCS(midRA, midDec)
        ax = WCSAxes(fig, [0.16, 0.1, 0.8, 0.8], wcs=wcs)
        fig.add_axes(ax)
    else:
        ax = plt.gca()

    # Plot facets
    markers = []
    for direction in directions_list:
        vertices = read_vertices(direction.vertices_file)
        RAverts = vertices[0]
        Decverts = vertices[1]
        xverts, yverts = factor.directions.radec2xy(RAverts, Decverts,
            refRA=midRA, refDec=midDec)
        xyverts = [np.array([xp, yp]) for xp, yp in zip(xverts, yverts)]
        mpl_poly = Polygon(np.array(xyverts), edgecolor='#a9a9a9', facecolor='none',
            clip_box=ax.bbox, picker=3.0, linewidth=2)
        mpl_poly.facet_name = direction.name
        mpl_poly.completed_ops = get_completed_ops(direction)
        mpl_poly.started_ops = get_started_ops(direction)
        if 'facetselfcal' in mpl_poly.completed_ops:
            if verify_subtract(direction):
                mpl_poly.set_edgecolor('g')
            else:
                mpl_poly.set_edgecolor('r')
        elif 'facetselfcal' in mpl_poly.started_ops:
            mpl_poly.set_edgecolor('y')
        elif 'facetimage' in mpl_poly.completed_ops:
            mpl_poly.set_edgecolor('b')
        elif 'facetimage' in mpl_poly.started_ops:
            mpl_poly.set_edgecolor('y')
        mpl_poly.selfcal_images = find_selfcal_images(direction)
        mpl_poly.facet_image = find_facet_image(direction)
        ax.add_patch(mpl_poly)

        # Show facet names
        xmid = (np.max(xverts) - np.min(xverts)) / 2.0 + np.min(xverts)
        ymid = (np.max(yverts) - np.min(yverts)) / 2.0 + np.min(yverts)
        marker = ax.text(xmid, ymid, direction.name, color='k', clip_on=True)
        markers.append(marker)

    ax.relim()
    ax.autoscale()
    ax.set_aspect('equal')

    if hasWCSaxes:
        RAAxis = ax.coords['ra']
        RAAxis.set_axislabel('RA', minpad=0.75)
        RAAxis.set_major_formatter('hh:mm:ss')
        DecAxis = ax.coords['dec']
        DecAxis.set_axislabel('Dec', minpad=0.75)
        DecAxis.set_major_formatter('dd:mm:ss')
        ax.coords.grid(color='black', alpha=0.5, linestyle='solid')
    else:
        plt.xlabel("RA (arb. units)")
        plt.ylabel("Dec (arb. units)")

    # Define coodinate formater to show RA and Dec under mouse pointer
    ax.format_coord = formatCoord

    # Show legend
    not_processed_patch = plt.Rectangle((0, 0), 1, 1, edgecolor='#a9a9a9',
        facecolor='none', linewidth=2)
    selfcal_ok_patch = plt.Rectangle((0, 0), 1, 1, edgecolor='g',
        facecolor='none', linewidth=2)
    selfcal_not_ok_patch =plt.Rectangle((0, 0), 1, 1, edgecolor='r',
        facecolor='none', linewidth=2)
    image_ok_patch = plt.Rectangle((0, 0), 1, 1, edgecolor='b',
        facecolor='none', linewidth=2)
    ax.legend([not_processed_patch, selfcal_ok_patch, selfcal_not_ok_patch, image_ok_patch],
              ['Unprocessed', 'Facetselfcal complete', 'Facetsefcal failed', 'Facetimage complete'])

    # Add check for mouse clicks
    fig.canvas.mpl_connect('pick_event', on_pick)

    # Show plot
    plt.show()
    plt.close(fig)

    # Clean up any temp pyrap images
    if os.path.exists('/tmp/tempimage'):
        shutil.rmtree('/tmp/tempimage')


def on_pick(event):
    facet = event.artist
    if event.mouseevent.button == 1:
        # Print info on left click
        log.info('Current state of reduction for {}:'.format(facet.facet_name))
        log.info('    Completed operations: {}'.format(facet.completed_ops))
        started_but_not_completed_ops = [op for op in facet.started_ops if
            not op in facet.completed_ops]
        log.info('      Running operations: {}'.format(started_but_not_completed_ops))

    if event.mouseevent.button == 3:
        # Open images (if any) on right click
        if os.path.exists('/tmp/tempimage'):
            shutil.rmtree('/tmp/tempimage')
        if len(facet.selfcal_images) > 0:
            log.info('Opening selfcal images for {}...'.format(facet.facet_name))
            im = pim.image(facet.selfcal_images)
            im.view()
        else:
            log.info('No selfcal images exist for {}'.format(facet.facet_name))
        if len(facet.facet_image) > 0:
            log.info('Opening facet image for {}...'.format(facet.facet_name))
            im2 = pim.image(facet.facet_image[0])
            im2.view()
        else:
            log.info('No full image of facet exists for {}'.format(facet.facet_name))

    plt.draw()


def get_completed_ops(direction):
    """
    Returns list of completed operations
    """
    has_state = direction.load_state()
    if has_state:
        return list(set(direction.completed_operations))
    else:
        return []


def get_started_ops(direction):
    """
    Returns list of started operations
    """
    has_state = direction.load_state()
    if has_state:
        return list(set(direction.started_operations))
    else:
        return []


def find_selfcal_images(direction):
    """
    Returns the filenames of selfcal images
    """
    selfcal_dir = os.path.join(direction.working_dir, 'results', 'facetselfcal',
        direction.name)
    if os.path.exists(selfcal_dir):
        selfcal_images = glob.glob(selfcal_dir+'/*.casa_image[0123]2.image.tt0')
        selfcal_images += glob.glob(selfcal_dir+'/*.casa_image42_iter*.image.tt0')
        if len(selfcal_images) == 0:
            selfcal_images = glob.glob(selfcal_dir+'/*.casa_image[0123]2.image')
            selfcal_images += glob.glob(selfcal_dir+'/*.casa_image42_iter*.image')
        selfcal_images.sort()
    else:
        selfcal_images = []

    return selfcal_images


def find_facet_image(direction):
    """
    Returns the filename of full facet image
    """
    selfcal_dir = os.path.join(direction.working_dir, 'results', 'facetselfcal',
        direction.name)
    image_dir = os.path.join(direction.working_dir, 'results', 'facetimage',
        direction.name)

    # Check selfcal and image directories. An image in the image directory is
    # preferred
    facet_image = []
    for d in [image_dir, selfcal_dir]:
        if os.path.exists(d) and len(facet_image) == 0:
            facet_image = glob.glob(d+'/*.wsclean_image_full2-MFS-image.fits')
            if len(facet_image) == 0:
                facet_image = glob.glob(d+'/*.casa_image_full2.image.tt0')
                if len(facet_image) == 0:
                    facet_image = glob.glob(d+'/*.casa_image_full2.image')

    return facet_image


def formatCoord(x, y):
    """Custom coordinate format"""
    global midRA, midDec
    RA, Dec = factor.directions.xy2radec([x], [y], midRA, midDec)
    return 'RA = {0:.2f} Dec = {1:.2f}'.format(RA[0], Dec[0])


def read_vertices(filename):
    """
    Returns facet vertices
    """
    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def verify_subtract(direction):
    """
    Checks selfcal success
    """
    verify_subtract_mapfile = os.path.join(direction.working_dir, 'results', 'facetselfcal',
        direction.name, 'mapfiles', 'verify_subtract.break.mapfile')
    if os.path.exists(verify_subtract_mapfile):
        ok_mapfile = DataMap.load(verify_subtract_mapfile)
        ok_flags = [ast.literal_eval(item.file) for item in ok_mapfile]
        if all(ok_flags):
            return True
        else:
            return False
    else:
        return False


def resample(array, factor):
    """
    Return resampled version of an image
    """

    nx, ny = np.shape(array)

    nx_new = nx // factor
    ny_new = ny // factor

    array2 = np.zeros((nx_new, ny))
    for i in range(nx_new):
        array2[i, :] = np.mean(array[i * factor:(i + 1) * factor, :], axis=0)

    array3 = np.zeros((nx_new, ny_new))
    for j in range(ny_new):
        array3[:, j] = np.mean(array2[:, j * factor:(j + 1) * factor], axis=1)

    return array3


def parset_read(parset_file):
    """
    Read and return subset of parset parameters

    Parameters
    ----------
    parset_file : str
        Filename of Factor-formated parset file

    Returns
    -------
    parset_dict : dict
        Dict of parset parameters

    """
    if not os.path.isfile(parset_file):
        log.error("Missing parset file (%s), I don't know what to do :'(" % (parset_file))
        sys.exit(1)

    log.info("Reading parset file: %s" % (parset_file))
    parset = ConfigParser.RawConfigParser()
    parset.read(parset_file)
    parset_dict = parset._sections['global']

    if not os.path.isdir(parset_dict['dir_working']):
        log.error('Cannot find the working directory. Please check the parset')
        sys.exit(1)

    if 'directions' in parset._sections.keys():
        parset_dict.update({'direction_specific': parset._sections['directions']})
    else:
        parset_dict.update({'direction_specific': {}})

    return parset_dict
