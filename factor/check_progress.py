"""
Module that checks and displays the progress of a run
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
from astropy.coordinates import Angle
from shapely.geometry import Polygon as SPolygon
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.parset import Parset
from factor.lib.direction import Direction
try:
    from matplotlib import pyplot as plt
    from matplotlib.patches import Polygon, Circle
    from matplotlib.ticker import FuncFormatter
    from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
except Exception as e:
    raise ImportError('PyPlot could not be imported. Plotting is not '
        'available: {0}'.format(e.message))
try:
    from wcsaxes import WCSAxes
    hasWCSaxes = True
except:
    hasWCSaxes = False

log = logging.getLogger('factor:progress')


def run(parset_file, trim_names=True):
    """
    Displays the current progress of a run
    """
    global all_directions

    log.info('Plotting facets...')
    log.info('Left-click on a facet to see its current state')
    log.info('Middle-click on a facet to display its image')
    log.info('Right-click on a facet to display its selfcal solutions and images')
    log.info('(In all cases, pan/zoom mode must be off)')
    log.info('Press "u" to update display (display is updated automatically every minute)')

    logging.root.setLevel(logging.ERROR)
    log.setLevel(logging.ERROR)
    all_directions = load_directions(parset_file)
    if len(all_directions) == 0:
        log.error('No directions found. Please check parset')
        sys.exit(1)
    plot_state(all_directions, trim_names=trim_names)


def load_directions(parset_file):
    """
    Return directions for a run
    """
    # Read parset
    parset = factor.parset.parset_read(parset_file, use_log_file=False)

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
    else:
        log.error('No directions found. Please run this tool after '
            'the directions have been defined')
        sys.exit(1)

    # Add the target to the directions list if desired
    target_ra = dir_parset['target_ra']
    target_dec = dir_parset['target_dec']
    target_radius_arcmin = dir_parset['target_radius_arcmin']
    target_has_own_facet = dir_parset['target_has_own_facet']
    if target_has_own_facet:
        if target_ra is not None and target_dec is not None and target_radius_arcmin is not None:
            # Make target object
            target = Direction('target', target_ra, target_dec,
                factor_working_dir=parset['dir_working'])

            # Add target to directions list
            directions.append(target)
        else:
            log.critical('target_has_own_facet = True, but target RA, Dec, or radius not found in parset')
            sys.exit(1)

    for direction in directions:
        has_state = direction.load_state()
        if has_state:
            direction_list.append(direction)

    return direction_list


def plot_state(directions_list, trim_names=True):
    """
    Plots the facets of a run
    """
    global midRA, midDec, fig, at

    # Set up coordinate system and figure
    points, midRA, midDec = factor.directions.getxy(directions_list)
    fig = plt.figure(1, figsize=(10,9))
    if hasWCSaxes:
        wcs = factor.directions.makeWCS(midRA, midDec)
        ax = WCSAxes(fig, [0.16, 0.1, 0.8, 0.8], wcs=wcs)
        fig.add_axes(ax)
    else:
        ax = plt.gca()

    field_x = min(points[0])
    field_y = max(points[1])
    adjust_xy = True
    while adjust_xy:
        adjust_xy = False
        for xy in points:
            dist = np.sqrt( (xy[0] - field_x)**2 + (xy[1] - field_y)**2 )
            if dist < 10.0:
                field_x -= 1
                field_y += 1
                adjust_xy = True
                break
    field_ra, field_dec = factor.directions.xy2radec([field_x], [field_y],
        refRA=midRA, refDec=midDec)
    field = Direction('field', field_ra[0], field_dec[0],
        factor_working_dir=directions_list[0].working_dir)
    directions_list.append(field)

    ax.set_title('Overview of FACTOR run in\n{}'.format(directions_list[0].working_dir))

    # Plot facets
    markers = []
    for direction in directions_list:
        if direction.name != 'field':
            vertices = read_vertices(direction.vertices_file)
            RAverts = vertices[0]
            Decverts = vertices[1]
            xverts, yverts = factor.directions.radec2xy(RAverts, Decverts,
                refRA=midRA, refDec=midDec)
            xyverts = [np.array([xp, yp]) for xp, yp in zip(xverts, yverts)]
            mpl_poly = Polygon(np.array(xyverts), edgecolor='#a9a9a9', facecolor='#F2F2F2',
                clip_box=ax.bbox, picker=3.0, linewidth=2)
        else:
            xverts = [field_x]
            yverts = [field_y]
            mpl_poly = Circle((field_x, field_y), radius=5.0, edgecolor='#a9a9a9', facecolor='#F2F2F2',
                clip_box=ax.bbox, picker=3.0, linewidth=2)
        mpl_poly.facet_name = direction.name
        mpl_poly.completed_ops = get_completed_ops(direction)
        mpl_poly.started_ops = get_started_ops(direction)
        mpl_poly.current_op = get_current_op(direction)
        set_patch_color(mpl_poly, direction)
        ax.add_patch(mpl_poly)

        # Add facet names
        if direction.name != 'field':
            poly_tuple = tuple([(xp, yp) for xp, yp in zip(xverts, yverts)])
            xmid = SPolygon(poly_tuple).centroid.x
            ymid = SPolygon(poly_tuple).centroid.y
        else:
            xmid = field_x
            ymid = field_y
        if trim_names:
            name = direction.name.split('_')[-1]
        else:
            name = direction.name
        marker = ax.text(xmid, ymid, name, color='k', clip_on=True,
            clip_box=ax.bbox, ha='center', va='bottom')
        markers.append(marker)

    # Add info box
    at = AnchoredText("Selected direction: None", prop=dict(size=12), frameon=True, loc=3)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)

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
        facecolor='#F2F2F2', linewidth=2)
    processing_patch = plt.Rectangle((0, 0), 1, 1, edgecolor='y',
        facecolor='#F2F5A9', linewidth=2)
    selfcal_ok_patch = plt.Rectangle((0, 0), 1, 1, edgecolor='g',
        facecolor='#A9F5A9', linewidth=2)
    selfcal_not_ok_patch =plt.Rectangle((0, 0), 1, 1, edgecolor='r',
        facecolor='#F5A9A9', linewidth=2)
    ax.legend([not_processed_patch, processing_patch, selfcal_ok_patch, selfcal_not_ok_patch],
              ['Unprocessed', 'Processing', 'Completed', 'Failed'])

    # Add check for mouse clicks and key presses
    fig.canvas.mpl_connect('pick_event', on_pick)
    fig.canvas.mpl_connect('key_press_event', on_press)

    # Add timer to update the plot every 60 seconds
    timer = fig.canvas.new_timer(interval=60000)
    timer.add_callback(update_plot)
    timer.start()

    # Show plot
    plt.show()
    plt.close(fig)

    # Clean up any temp pyrap images
    if os.path.exists('/tmp/tempimage'):
        shutil.rmtree('/tmp/tempimage')


def on_pick(event):
    """
    Handle picks with the mouse
    """
    global all_directions, at, fig

    facet = event.artist
    direction = None
    if hasattr(facet, 'facet_name'):
        for d in all_directions:
            if d.name == facet.facet_name:
                direction = d
                break
        if direction is None:
            return
    else:
        return

    if event.mouseevent.button == 1: # left click
        # Print info
        current_op = get_current_op(direction)
        info = 'Selected direction: {}\n'.format(facet.facet_name)
        completed_ops = get_completed_ops(direction)
        if len(completed_ops) == 0:
            info += 'Completed ops: None\n'
        else:
            info += 'Completed ops: {}\n'.format(', '.join(completed_ops))
        info += 'Current op: {}'.format(current_op)
        if current_op is not None:
            current_step, current_index, num_steps, start_time = get_current_step(direction)
            if current_step is not None:
                info += '\n- Started at: {}\n'.format(start_time)
                info += '- Current step: {0} (step {1} of {2})'.format(
                    current_step, current_index+1, num_steps)

    elif event.mouseevent.button == 2: # middle click
        # Open full facet image (if any)
        if os.path.exists('/tmp/tempimage'):
            shutil.rmtree('/tmp/tempimage')
        facet_image = find_facet_image(direction)
        if len(facet_image) > 0:
            info = 'Opening facet image for {}...'.format(facet.facet_name)
            im2 = pim.image(facet_image[0])
            im2.view()
        else:
            info = 'No image of facet exists for {}'.format(facet.facet_name)

    elif event.mouseevent.button == 3: # right click
        # Open selfcal images (if any)
        if os.path.exists('/tmp/tempimage'):
            shutil.rmtree('/tmp/tempimage')
        selfcal_images = find_selfcal_images(direction)
        if len(selfcal_images) > 0:
            info = 'Opening selfcal images for {}...'.format(facet.facet_name)
            im = pim.image(selfcal_images)
            im.view()
        else:
            info = 'No selfcal images exist for {}'.format(facet.facet_name)

        # Open parmdbplot of selfcal instrument table (if any)
        selfcal_parmdb = find_selfcal_parmdb(direction)
        if selfcal_parmdb is not None:
            info += '\nOpening final selfcal solutions for {}...'.format(facet.facet_name)
            os.system('parmdbplot.py {} &'.format(selfcal_parmdb))
        else:
            info += '\nFinal selfcal solutions do not exist for {}'.format(facet.facet_name)

    else:
        return

    # Update text box
    c = at.get_child()
    c.set_text(info)

    # Update patch color
    set_patch_color(facet, direction)

    # Redraw
    fig.canvas.draw()


def on_press(event):
    """
    Handle key presses
    """
    global fig, all_directions, at

    if event.key == 'u':
        info = 'Updating display...'
        c = at.get_child()
        c.set_text(info)
        fig.canvas.draw()
        update_plot()
        info += '\n...done'
        c = at.get_child()
        c.set_text(info)
        fig.canvas.draw()


def update_plot():
    """
    Update the plot
    """
    global fig, all_directions, at

    ax = plt.gca()
    for a in ax.patches:
        if hasattr(a, 'facet_name'):
            for d in all_directions:
                if d.name == a.facet_name:
                    set_patch_color(a, d)
                    a.selfcal_images = find_selfcal_images(d)
                    a.facet_image = find_facet_image(d)
    fig.canvas.draw()


def set_patch_color(a, d):
    """
    Sets face and edge color of patch
    """
    a.completed_ops = get_completed_ops(d)
    a.started_ops = get_started_ops(d)
    a.current_op = get_current_op(d)
    if a.current_op is not None:
        a.set_edgecolor('y')
        a.set_facecolor('#F2F5A9')
    elif 'facetselfcal' in a.completed_ops:
        if verify_subtract(d):
            a.set_edgecolor('g')
            a.set_facecolor('#A9F5A9')
        else:
            a.set_edgecolor('r')
            a.set_facecolor('#F5A9A9')
    elif len(a.completed_ops) > 0:
        a.set_edgecolor('g')
        a.set_facecolor('#A9F5A9')


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


def get_current_op(direction):
    """
    Returns name of current op (assumes there is just one)
    """
    started_ops = get_started_ops(direction)
    completed_ops = get_completed_ops(direction)
    started_but_not_completed_ops = [op for op in started_ops if not op in completed_ops]

    if len(started_but_not_completed_ops) == 0:
        return None
    else:
        return started_but_not_completed_ops[0]


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


def find_selfcal_parmdb(direction):
    """
    Returns the filename of selfcal parmdb
    """
    selfcal_dir = os.path.join(direction.working_dir, 'results', 'facetselfcal',
        direction.name)
    if os.path.exists(selfcal_dir):
        selfcal_parmdb = glob.glob(selfcal_dir+'/*.merge_selfcal_parmdbs')
        if len(selfcal_parmdb) == 0:
            selfcal_parmdb = None
        else:
            selfcal_parmdb = selfcal_parmdb[0]
    else:
        selfcal_parmdb = None

    return selfcal_parmdb


def find_facet_image(direction):
    """
    Returns the filename of full facet image
    """
    selfcal_dir = os.path.join(direction.working_dir, 'results', 'facetselfcal',
        direction.name)
    image_dir = os.path.join(direction.working_dir, 'results', 'facetimage',
        direction.name)

    # Check selfcal and image directories. An image in the facetimage directory
    # is preferred
    facet_image = []
    for d in [image_dir, selfcal_dir]:
        if os.path.exists(d) and len(facet_image) == 0:
            facet_image = glob.glob(d+'/*.wsclean_image_full2-MFS-image.fits')
            if len(facet_image) == 0:
                facet_image = glob.glob(d+'/*.wsclean_image_full2-image.fits')
                if len(facet_image) == 0:
                    facet_image = glob.glob(d+'/*.casa_image_full2.image.tt0')
                    if len(facet_image) == 0:
                        facet_image = glob.glob(d+'/*.casa_image_full2.image')

    return facet_image


def get_current_step(direction):
    """
    Returns name and index of current step, total number of steps, and start time
    of current operation
    """
    current_op = get_current_op(direction)
    if current_op is None:
        return (None, None, None, None)

    statefile = os.path.join(direction.working_dir, 'results', current_op,
        direction.name, 'statefile')
    if os.path.exists(statefile):
        f = open(statefile, 'r')
        d = pickle.load(f)
        f.close()
    else:
        return (None, None, None, None)

    mapfiles = set([s[1]['mapfile'] for s in d[1]])
    current_index = len(mapfiles)
    current_steps = get_current_op_step_names(direction)
    if current_index >= len(current_steps):
        current_index = len(current_steps) - 1
    start_time = d[0]['start_time']

    return (current_steps[current_index], current_index, len(current_steps), start_time)


def get_current_op_step_names(direction):
    """
    Returns lits of step names for current operation
    """
    current_op = get_current_op(direction)
    parset_file = os.path.join(direction.working_dir, 'results', current_op,
        direction.name, 'pipeline.parset')
    parset = Parset()
    parset.adoptFile(parset_file)
    pipeline_args = parset.makeSubset(parset.fullModuleName('pipeline') + '.')
    step_name_list = pipeline_args.getStringVector('steps')

    # Filter out plugin steps
    filter_step_name_list = []
    for stepname in step_name_list:
        fullparset = parset.makeSubset(parset.fullModuleName(str(stepname)) + '.')
        subparset = fullparset.makeSubset(fullparset.fullModuleName('control') + '.')
        try:
            kind_of_step = subparset.getString('kind')
        except:
            kind_of_step = 'recipe'
        if kind_of_step != 'plugin':
            if kind_of_step == 'loop':
                loopsteps = subparset.getStringVector('loopsteps')
                for loopstep in loopsteps:
                    fullparset_loop = parset.makeSubset(parset.fullModuleName(str(loopstep)) + '.')
                    subparset_loop = fullparset_loop.makeSubset(fullparset_loop.fullModuleName('control') + '.')
                    try:
                        kind_of_loop_step = subparset_loop.getString('kind')
                    except:
                        kind_of_loop_step = 'recipe'
                    if kind_of_loop_step != 'plugin':
                        filter_step_name_list.append(loopstep)
            else:
                filter_step_name_list.append(stepname)

    return filter_step_name_list


def formatCoord(x, y):
    """Custom coordinate format"""
    global midRA, midDec

    RA, Dec = factor.directions.xy2radec([x], [y], midRA, midDec)
    RA_str = Angle(RA[0], unit='deg').to_string('hour')
    Dec_str = Angle(Dec[0], unit='deg').to_string('deg')
    return 'RA = {0} Dec = {1}'.format(RA_str, Dec_str)


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
