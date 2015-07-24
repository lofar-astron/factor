"""
Module that holds all direction-related functions
"""
import os
import numpy as np
import logging
from factor.lib.direction import Direction
import sys
from scipy.spatial import Delaunay


log = logging.getLogger('factor.directions')


def directions_read(directions_file, factor_working_dir):
    """
    Read a Factor-formatted directions file and return list of Direction objects

    Parameters
    ----------
    directions_file : str
        Filename of Factor-formated directions file
    factor_working_dir : str
        Full path of working directory

    Returns
    -------
    directions : list of Direction objects
        List of Direction objects

    """
    from astropy.coordinates import Angle

    if not os.path.isfile(directions_file):
        log.critical("Directions file (%s) not found." % (directions_file))
        sys.exit(1)

    log.info("Reading directions file: %s" % (directions_file))
    try:
        types = np.dtype({'names': ['name', 'radec', 'atrous_do', 'mscale_field_do',
            'cal_imsize', 'solint_p', 'solint_a', 'field_imsize', 'dynamic_range',
            'region_selfcal', 'region_field', 'peel_skymodel', 'outlier_source',
            'cal_size_deg', 'cal_flux_jy'], 'formats':['S255', 'S255', 'S5',
            'S5', int, int, int, int, 'S2', 'S255', 'S255', 'S255',
            'S5', float, float]})
        directions = np.genfromtxt(directions_file, comments='#', dtype=types)
    except ValueError:
        types = np.dtype({'names': ['name', 'radec', 'atrous_do', 'mscale_field_do',
            'cal_imsize', 'solint_p', 'solint_a', 'field_imsize', 'dynamic_range',
            'region_selfcal', 'region_field', 'peel_skymodel', 'outlier_source'],
            'formats':['S255', 'S255', 'S5',
            'S5', int, int, int, int, 'S2', 'S255', 'S255', 'S255',
            'S5', float, float]})
        directions = np.genfromtxt(directions_file, comments='#', dtype=types)

    data = []
    for direction in directions:
        RAstr, Decstr = direction['radec'].split(',')
        ra = Angle(RAstr).to('deg').value
        dec = Angle(Decstr).to('deg').value

        # some checks on values
        if np.isnan(ra) or ra < 0 or ra > 360:
            log.error('RA %f is wrong for direction: %s. Ignoring direction.'
            	% (direction['radec'], direction['name']))
            continue
        if np.isnan(dec) or dec < -90 or dec > 90:
            log.error('DEC %f is wrong for direction: %s. Ignoring direction.'
                % (direction['radec'], direction['name']))
            continue
        if direction['atrous_do'].lower() == 'true':
            atrous_do = True
        else:
            atrous_do = False
        if direction['mscale_field_do'].lower() == 'true':
            mscale_field_do = True
        else:
            mscale_field_do = False
        if direction['outlier_source'].lower() == 'true':
            outlier_source = True
        else:
            outlier_source = False
        if (direction['solint_a'] <= 0 or direction['solint_p'] <= 0) and \
            np.isnan(direction['apparent_flux']):
            log.error('One of more of the solution intervals is invalid and no '
                'apparent flux is specified for direction {0}. Ignoring '
                'direction.'.format(direction['name']))
            continue

        # set defaults
        if direction['solint_a'] <= 0:
            direction['solint_a'] = 60
        if direction['solint_p'] <= 0:
            direction['solint_p'] = 1
        if len(direction) > 13:
            if direction['cal_size_deg'] <= 0.0 or np.isnan(direction['cal_size_deg']):
                cal_size_deg = None
            else:
                cal_size_deg = direction['cal_size_deg']
            if np.isnan(direction['cal_flux_jy']):
                cal_flux_jy = None
            else:
                cal_flux_jy = direction['cal_flux_jy']
        else:
            cal_size_deg = None
            cal_flux_jy = None

        data.append(Direction(direction['name'], ra, dec,
        	atrous_do, mscale_field_do,
        	direction['cal_imsize'], direction['solint_p'],
        	direction['solint_a'], direction['field_imsize'],
        	direction['dynamic_range'], direction['region_selfcal'],
        	direction['region_field'], direction['peel_skymodel'],
        	outlier_source, factor_working_dir, False,
        	cal_size_deg, cal_flux_jy))

    return data


def make_directions_file_from_skymodel(bands, flux_min_Jy, size_max_arcmin,
    directions_separation_max_arcmin, directions_max_num=None,
    interactive=False):
    """
    Selects appropriate calibrators from sky models and makes the directions file

    Parameters
    ----------
    bands : list of Band objects
        Input Bands with dir-indep skymodels
    flux_min_Jy : float
        Minimum flux for a calibrator in Jy
    size_max_arcmin : float
        Maximum size for a calibrator in arcmin
    directions_separation_max_arcmin : float
        Maximum separation in arcmin between two calibrators for gouping into a
        single direction
    directions_max_num : int, optional
        Limit total number of directions to this value
    interactive : bool, optional
        If True, plot the directions and ask for approval

    Returns
    -------
    directions_file : str
        Filename of resulting Factor-formatted directions file

    """
    # Use sky model of lowest-frequency band
    freqs = [band.freq for band in bands]
    min_freq_indx = np.argmin(freqs)
    band = bands[min_freq_indx]
    directions_file = 'factor_directions.txt'
    ds9_directions_file = 'factor_directions_ds9.reg'

    # Load sky model and filter it
    s = make_initial_skymodel(band)

    # Filter larger patches
    sizes = s.getPatchSizes(units='arcmin', weight=True)
    s.select(sizes < size_max_arcmin, aggregate=True, force=True)
    if len(s) == 0:
        log.critical("No sources found that meet the specified max size criteria.")
        sys.exit(1)
    log.info('Found {0} sources with sizes below {1} '
        'arcmin'.format(len(s.getPatchNames()), size_max_arcmin))

    # Look for nearby pairs
    log.info('Merging sources within {0} arcmin of each other...'.format(
        directions_separation_max_arcmin))
    pRA, pDec = s.getPatchPositions(asArray=True)
    for ra, dec in zip(pRA.tolist()[:], pDec.tolist()[:]):
        dist = s.getDistance(ra, dec, byPatch=True, units='arcmin')
        nearby = np.where(dist < directions_separation_max_arcmin)
        if len(nearby[0]) > 1:
            patches = s.getPatchNames()[nearby]
            s.merge(patches.tolist())

    # Filter fainter patches
    s.select('I > {0} Jy'.format(flux_min_Jy), aggregate='sum', force=True)
    if len(s) == 0:
        log.critical("No sources found that meet the specified min flux criteria.")
        sys.exit(1)
    log.info('Found {0} sources with fluxes above {1} Jy'.format(
        len(s.getPatchNames()), flux_min_Jy))

    # Trim directions list to get directions_total_num of directions
    if directions_max_num is not None:
        dir_fluxes = s.getColValues('I', aggregate='sum')
        dir_fluxes_sorted = dir_fluxes.tolist()
        dir_fluxes_sorted.sort(reverse=True)
        cut_jy = dir_fluxes_sorted[-1]
        while len(dir_fluxes_sorted) > directions_max_num:
            cut_jy = dir_fluxes_sorted.pop() + 0.00001
        s.remove('I < {0} Jy'.format(cut_jy), aggregate='sum')

    log.info('Kept {0} directions in total'.format(len(s.getPatchNames())))

    # Write the file
    s.setPatchPositions(method='mid')
    log.info("Writing directions file: %s" % (directions_file))
    s.write(fileName=directions_file, format='factor', sortBy='I', clobber=True)

    return directions_file


def make_initial_skymodel(band):
    """
    Makes the initial skymodel used to adjust facet edges

    Parameters
    ----------
    band : Band object
        Band to use for sky model generation

    Returns
    -------
    s : LSMTool Skymodel object
        Resulting sky model

    """
    import lsmtool

    # Load sky model and filter it
    s = lsmtool.load(band.skymodel_dirindep)

    # Group to clean components by thresholding after convolving model with
    # 1 arcmin beam
    s.group('threshold', FWHM='60.0 arcsec', root='facet')
    s.remove('Patch = patch_*', force=True) # Remove sources that did not threshold
    if len(s) == 0:
        log.critical("No sources found through thresholding.")
        sys.exit(1)
    log.info('Found {0} sources through thresholding'.format(
        len(s.getPatchNames())))

    # Filter out sources that lie outside of FWHM of FOV
    if not hasattr(band, 'fwhm_deg'):
        band.set_image_sizes()
    log.info('Removing sources beyond 2 * FWHM of the primary beam...')
    dist = s.getDistance(band.ra, band.dec, byPatch=True)
    s.remove(dist > band.fwhm_deg, aggregate=True)

    # Save this sky model for later checks of sources falling on facet edges
    s.write(fileName='results/initial.skymodel', clobber=True)

    return s


def group_directions(directions, one_at_a_time=True, n_per_grouping={'1':0},
    allow_reordering=True):
    """
    Sorts directions into groups that can be selfcaled simultaneously

    Directions are grouped by flux and then optionally reordered to maximize
    the miniumum separation between sources in a group

    Parameters
    ----------
    directions : list of Direction objects
        List of input directions to group
    one_at_a_time : bool, optional
        If True, run one direction at a time
    n_per_grouping : dict, optional
        Dict specifying the total number of sources at each grouping level. The
        sources at each grouping level can be reordered to maximize the
        minimum separation between sources within a group
    allow_reordering : bool, optional
        If True, allow sources in neighboring groups to be reordered to increase
        the minimum separation between sources within a group

    Returns
    -------
    direction_groups : list of lists
        List of direction groups

    """
    from random import shuffle

    direction_groups = []
    if one_at_a_time or n_per_grouping == {'1': 0}:
        for d in directions:
            direction_groups.append([d])
        log.debug('Processing each direction in series')
    else:
        def find_min_separation(group):
            """
            Finds the minimum separation in degrees between sources in a group
            """
            sep = []
            for direction1 in group:
                for direction2 in group:
                    if direction1 != direction2:
                        sep.append(calculateSeparation(direction1.ra, direction1.dec,
                            direction2.ra, direction2.dec))
            return min(sep)

        # Divide based on flux (assuming order is decreasing flux)
        log.info('Dividing directions into groups...')
        grouping_levels = [int(g) for g in n_per_grouping.iterkeys()]
        grouping_levels.sort()
        for i, g in enumerate(grouping_levels):
            if i == 0:
                start = 0
            else:
                start = end
            if i < len(grouping_levels)-1:
                if n_per_grouping[str(g)] <= 0:
                    # check for 0 as indicator of "take all the rest"
                    end = len(directions)
                else:
                    end = start + n_per_grouping[str(g)]
            else:
                end = len(directions)
            if end > len(directions):
                end = len(directions)
            if end > start:
                for j in range(start, end, g):
                    gstart = j
                    gend = j + g
                    if gend > end:
                        gend = end
                    if j == 0:
                        direction_groups = [directions[gstart: gend]]

                    else:
                        direction_groups += [directions[gstart: gend]]

        # Reorganize groups in each grouping level based on distance. This
        # is done by swapping the directions of neighboring groups randomly
        # and picking the group with the largest minimum separation
        if allow_reordering:
            log.info('Reordering directions to obtain max separation...')
            direction_groups_orig = direction_groups[:]
            if len(direction_groups_orig) > 1:
                for i in range(0, len(direction_groups_orig), 2):
                    group1 = direction_groups_orig[i]
                    if len(group1) > 1:
                        if i < len(direction_groups)-1:
                            k = i + 1
                        else:
                            k = i - 1
                        group2 = direction_groups_orig[k]

                        min_sep_global = 0.0 # degrees
                        for j in range(10):
                            group_merged = group1[:] + group2[:]
                            shuffle(group_merged)
                            group1_test = group_merged[0: len(group1)]
                            group2_test = group_merged[len(group1):]
                            min_sep1 = find_min_separation(group1_test)
                            min_sep2 = find_min_separation(group2_test)
                            min_sep = min(min_sep1, min_sep2)
                            if min_sep > min_sep_global:
                                min_sep_global = min_sep
                                group1_best = group1_test
                                group2_best = group2_test
                        direction_groups[i] = group1_best
                        direction_groups[k] = group2_best
        log.debug('Processing directions in the following groups:')
        for i, group in enumerate(direction_groups):
            log.debug('Group {0}: {1}'.format(i+1, [d.name for d in group]))

    return direction_groups


def thiessen(directions_list, bounds_scale=0.52, band=None, check_edges=False,
    target_ra=None, target_dec=None, target_radius_arcmin=None):
    """
    Return list of thiessen polygons and their widths in degrees

    Parameters
    ----------
    directions_list : list of Direction objects
        List of input directions
    bounds_scale : int, optional
        Scale to use for bounding box
    band : Band object, optional
        Band to use to check for source near facet edges
    check_edges : bool, optional
        If True, check whether any know source falls on a facet edge. If sources
        are found that do, the facet is adjusted
    target_ra : str, optional
        RA of target source. E.g., '14h41m01.884'
    target_dec : str, optional
        Dec of target source. E.g., '+35d30m31.52'
    target_radius_arcmin : float, optional
        Radius in arcmin of target source

    Returns
    -------
    thiessen_polys_deg : list
        List of polygon RA and Dec vertices in degrees
    width_deg : list
        List of polygon bounding box width in degrees

    """
    import lsmtool
    try:
        import shapely.geometry
        from shapely.ops import cascaded_union
        has_shapely = True
    except ImportError:
        log.warn('Shapely could not be imported. Facet polygons will not be '
            'adjusted to avoid known sources.')
        has_shapely = False
    from itertools import combinations
    from astropy.coordinates import Angle

    points, midRA, midDec = getxy(directions_list)
    points = points.T

    x_scale, y_scale = (points.min(axis=0) - points.max(axis=0)) * bounds_scale

    means = np.ones((32, 2)) * points.mean(axis=0)

    radius = np.sqrt(x_scale**2 + y_scale**2)
    angles = [np.pi/16.0*i for i in range(0, 32)]
    offsets = []
    for ang in angles:
        offsets.append([np.cos(ang), np.sin(ang)])
    scale_offsets = radius * np.array(offsets)
    outer_box = means + scale_offsets

    points = np.vstack([points, outer_box])
    tri = Delaunay(points)
    circumcenters = np.array([_circumcenter(tri.points[t])
                              for t in tri.vertices])
    thiessen_polys = [_thiessen_poly(tri, circumcenters, n)
                      for n in range(len(points) - 32)]

    # Check for sources near / on facet edges and adjust regions accordingly
    if has_shapely and check_edges:
        log.info('Adjusting facets to avoid sources...')
        if os.path.exists('results/initial.skymodel'):
            s = lsmtool.load('results/initial.skymodel')
        elif band is not None:
            s = make_initial_skymodel(band)
        else:
            log.error('A band must be given for edge checking')
            sys.exit(1)
        RA, Dec = s.getPatchPositions(asArray=True)
        sx, sy = radec2xy(RA, Dec, refRA=midRA, refDec=midDec)
        sizes = s.getPatchSizes(units='degree').tolist()

        if target_ra is not None and target_dec is not None and target_radius_arcmin is not None:
            log.info('Including target ({0}, {1}) in facet adjustment'.format(
                target_ra, target_dec))
            tra = Angle(target_ra).to('deg').value
            tdec = Angle(target_dec).to('deg').value
            tx, ty = radec2xy([tra], [tdec], refRA=midRA, refDec=midDec)
            sx.extend(tx)
            sy.extend(ty)
            sizes.append(target_radius_arcmin*2.0/1.2/60.0)

        # Filter sources to get only those close to a boundary. We need to iterate
        # until no sources are found
        niter = 0
        while niter < 3:
            niter += 1
            ind_near_edge = []
            for i, thiessen_poly in enumerate(thiessen_polys):
                polyv = np.vstack(thiessen_poly)
                poly_tuple = tuple([(x, y) for x, y in zip(polyv[:, 0], polyv[:, 1])])
                poly = Polygon(polyv[:, 0], polyv[:, 1])
                dists = poly.is_inside(sx, sy)
                for j, dist in enumerate(dists):
                    pix_radius = sizes[j] * 1.2 / 2.0 / 0.066667 # radius of source in pixels
                    if abs(dist) < pix_radius and j not in ind_near_edge:
                        ind_near_edge.append(j)
            if len(ind_near_edge) == 0:
                break
            sx_filt = np.array(sx)[ind_near_edge]
            sy_filt = np.array(sy)[ind_near_edge]
            sizes_filt = np.array(sizes)[ind_near_edge]

            # Adjust all facets for each source near a boundary
            for x, y, size in zip(sx_filt, sy_filt, sizes_filt):
                for i, thiessen_poly in enumerate(thiessen_polys):
                    polyv = np.vstack(thiessen_poly)
                    poly_tuple = tuple([(xp, yp) for xp, yp in zip(polyv[:, 0], polyv[:, 1])])
                    poly = Polygon(polyv[:, 0], polyv[:, 1])
                    dist = poly.is_inside(x, y)
                    p1 = shapely.geometry.Polygon(poly_tuple)

                    pix_radius = size * 1.2 / 2.0 / 0.066667 # size of source in pixels
                    if abs(dist) < pix_radius:
                        p2 = shapely.geometry.Point((x, y))
                        p2buf = p2.buffer(pix_radius)
                        if dist < 0.0:
                            # If point is outside, difference the polys
                            p1 = p1.difference(p2buf)
                        else:
                            # If point is inside, union the polys
                            p1 = p1.union(p2buf)
                        try:
                            xyverts = [np.array([xp, yp]) for xp, yp in
                                zip(p1.exterior.coords.xy[0].tolist(),
                                p1.exterior.coords.xy[1].tolist())]
                        except AttributeError:
                            log.error('Source avoidance has caused a facet to be '
                                'divided into multple parts. Please adjust the '
                                'parameters (e.g., if a target source is specified, '
                                'reduce its radius if possible)')
                            sys.exit(1)
                        thiessen_polys[i] = xyverts

    # Convert from x, y to RA, Dec and find width of facet and facet center
    for d, poly in zip(directions_list, thiessen_polys):
        poly = np.vstack([poly, poly[0]])
        ra, dec = xy2radec(poly[:, 0], poly[:, 1], midRA, midDec)
        thiessen_poly_deg = [np.array(ra[0: -1]), np.array(dec[0: -1])]

        # Find size and centers of regions in degrees
        xmin = np.min(poly[:, 0])
        xmax = np.max(poly[:, 0])
        xmid = xmin + int((xmax - xmin) / 2.0)
        ymin = np.min(poly[:, 1])
        ymax = np.max(poly[:, 1])
        ymid = ymin + int((ymax - ymin) / 2.0)

        ra1, dec1 = xy2radec([xmin], [ymin], midRA, midDec)
        ra2, dec2 = xy2radec([xmax], [ymax], midRA, midDec)
        ra3, dec3 = xy2radec([xmax], [ymin], midRA, midDec)
        ra_center, dec_center = xy2radec([xmid], [ymid], midRA, midDec)
        ra_width_deg = calculateSeparation(ra1, dec1, ra3, dec3)
        dec_width_deg = calculateSeparation(ra3, dec3, ra2, dec2)
        width_deg = max(ra_width_deg.value, dec_width_deg.value)

        d.vertices = thiessen_poly_deg
        d.width = width_deg
        d.facet_ra = ra_center[0]
        d.facet_dec = dec_center[0]


def make_region_file(vertices, outputfile):
    """
    Make a CASA region file for given vertices

    Parameters
    ----------
    vertices : list
        List of direction RA and Dec vertices in degrees
    outputfile : str
        Name of output region file

    Returns
    -------
    region_filename : str
        Name of region file

    """
    lines = ['#CRTFv0\n\n']
    xylist = []
    RAs = vertices[0][0:-1] # trim last point, as it is a repeat of the first
    Decs = vertices[1][0:-1]
    for x, y in zip(RAs, Decs):
        xylist.append('[{0}deg, {1}deg]'.format(x, y))
    lines.append('poly[{0}]\n'.format(', '.join(xylist)))

    with open(outputfile, 'wb') as f:
        f.writelines(lines)


def make_ds9_region_file(directions, outputfile):
    """
    Make a ds9 region file for given vertices and centers

    Parameters
    ----------
    directions : list
        List of Direction objects
    outputfile : str
        Name of output region file

    """
    lines = []
    lines.append('# Region file format: DS9 version 4.0\nglobal color=green '
                 'font="helvetica 10 normal" select=1 highlite=1 edit=1 '
                 'move=1 delete=1 include=1 fixed=0 source=1\nfk5\n')

    for direction in directions:
        xylist = []
        RAs = direction.vertices[0]
        Decs = direction.vertices[1]
        for x, y in zip(RAs, Decs):
            xylist.append('{0}, {1}'.format(x, y))
        lines.append('polygon({0})\n'.format(', '.join(xylist)))
        lines.append('point({0}, {1}) # point=cross width=2 text={{{2}}}\n'.
            format(direction.ra, direction.dec, direction.name))

    with open(outputfile, 'wb') as f:
        f.writelines(lines)


def make_ds9_calimage_file(directions, outputfile):
    """
    Make a ds9 image region file for given calibrator size

    Parameters
    ----------
    directions : list
        List of Direction objects
    outputfile : str
        Name of output region file

    """
    lines = []
    lines.append('# Region file format: DS9 version 4.0\nglobal color=yellow '
                 'font="helvetica 10 normal" select=1 highlite=1 edit=1 '
                 'move=1 delete=1 include=1 fixed=0 source=1\nfk5\n')

    for direction in directions:
        imsize = direction.cal_imsize * direction.cellsize_selfcal_deg * 3600 # arcsec
        imsize_unmasked = 0.8 * imsize
        RAs = direction.vertices[0]
        Decs = direction.vertices[1]
        lines.append('box({0}, {1}, {2}", {2}") # text={{{3}}}\n'.
            format(direction.ra, direction.dec, imsize, direction.name))
        lines.append('box({0}, {1}, {2}", {2}")\n'.
            format(direction.ra, direction.dec, imsize_unmasked))

    with open(outputfile, 'wb') as f:
        f.writelines(lines)


def plot_thiessen(directions_list, bounds_scale=2):
    """
    Plot thiessen polygons for a given set of points

    Parameters
    ----------
    directions_list : list of Direction objects
        List of input directions
    bounds_scale : int, optional
        Scale to use for bounding box

    """
    from matplotlib import pyplot as plt

    points, midRA, midDec = getxy(directions_list)
    points = points.T
    polys, _ = thiessen(directions_list, bounds_scale)
    plt.scatter(points[:, 0], points[:, 1])
    for poly in polys:
        poly = np.vstack([poly, poly[0]])
        plt.plot(poly[:, 0], poly[:, 1], 'r')
        poly = np.vstack([poly, poly[0]])
    plt.show()


def _any_equal(arr, n):
    """for a given Mx3 array, returns a 1xM array containing indices
    of rows where any of the columns are equal to n.
    """
    return np.where((arr[:, 0] == n) | (arr[:, 1] == n) | (arr[:, 2] == n))[0]


def _circumcenter(vertices):
    """returns the circumcenter of a triangle.
    ``vertices`` should be a np.array of size (3,2) containing the
    points of the triangle
    """
    ax, ay, bx, by, cx, cy = vertices.flatten()

    D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))

    # don't divide by 0
    if D == 0:
        D = 0.000000001

    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) +
        (cx**2 + cy**2) * (ay - by)) / D
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) +
        (cx**2 + cy**2) * (bx - ax)) / D

    return ux, uy


def _find_triangles_for_vertex(tri, n):
    """returns all of the indices of the triangles for the nth vertex
    of a given a scipy.spatial.Delaunay object
    """
    # grab the list of triangles that touch this vertex
    triangles = tri.vertices[_any_equal(tri.vertices, n)]

    # we want to sort the triangles so that neighbors are together,
    # just start with the first triangle
    sorted_triangles = [triangles[0]]

    # initialize values
    if triangles[0][0] != n:
        previous_vertex_idx = triangles[0][0]
    else:
        previous_vertex_idx = triangles[0][1]

    # just stash the common vertex for checking if we're sorted
    # clockwise later on
    common_edge_vertex_idx = previous_vertex_idx

    # loop through the triangles; previous_vertex_index will be the
    # index to the vertex we used in the previous triangle
    for i in triangles[1:]:
        this_triangle = sorted_triangles[-1]

        # find the vertex of the triangle that is not the central
        # vertex and is not shared with the previous triangle
        next_vertex_idx = this_triangle[(this_triangle != n) & (this_triangle
            != previous_vertex_idx)]

        # append the next triangle (note: match will return both the
        # previous triangle and the next triangle, since they both
        # contain the shared vertex)
        matching_triangles = triangles[_any_equal(triangles, next_vertex_idx)]
        if np.all(this_triangle == matching_triangles[0]):
            sorted_triangles.append(matching_triangles[1])
        else:
            sorted_triangles.append(matching_triangles[0])

        previous_vertex_idx = next_vertex_idx

    sorted_triangle_indices = [
        int(np.where(np.all(tri.vertices[:] == triangle, axis=1))[0])
        for triangle in sorted_triangles]

    # if we're sorted counter-clockwise, then we need to reverse order
    test_point = tri.points[triangles[0][(triangles[0] != n) & (triangles[0]
        != common_edge_vertex_idx)]].flatten()
    if not _is_right(tri.points[n], tri.points[common_edge_vertex_idx], test_point):
        return sorted_triangle_indices[::-1]

    # otherwise we're good
    return sorted_triangle_indices


def _is_right(a, b, p):
    """given a line (defined by points a and b) and a point (p),
    return true if p is to the right of the line and false otherwise
    raises a ValueError if p lies is colinear with a and b
    """
    ax, ay = a[0], a[1]
    bx, by = b[0], b[1]
    px, py = p[0], p[1]
    value = (bx - ax) * (py - ay) - (by - ay) * (px - ax)

    if value == 0:
        raise ValueError(
            "p is colinear with a and b, 'tis neither right nor left.")

    return value < 0


def _thiessen_poly(tri, circumcenters, n):
    """given a Delaunay triangulation object, calculates a thiessen
    polygon for the vertex index n
    """
    triangles = _find_triangles_for_vertex(tri, n)
    triangles = np.hstack((triangles, triangles[0]))
    return [circumcenters[t] for t in triangles]


def getxy(directions_list):
    """
    Returns array of projected x and y values.

    Parameters
    ----------
    directions_list : list
        List of direction objects

    Returns
    -------
    x, y, midRA, midDec : numpy array, numpy array, float, float
        arrays of x and y values and the midpoint RA and
        Dec values

    """
    import numpy as np

    if len(directions_list) == 0:
        return np.array([0, 0]), 0, 0

    RA = []
    Dec = []
    for direction in directions_list:
        RA.append(direction.ra)
        Dec.append(direction.dec)
    x, y  = radec2xy(RA, Dec)

    # Refine x and y using midpoint
    if len(x) > 1:
        xmid = min(x) + (max(x) - min(x)) / 2.0
        ymid = min(y) + (max(y) - min(y)) / 2.0
        xind = np.argsort(x)
        yind = np.argsort(y)
        try:
            midxind = np.where(np.array(x)[xind] > xmid)[0][0]
            midyind = np.where(np.array(y)[yind] > ymid)[0][0]
            midRA = RA[xind[midxind]]
            midDec = Dec[yind[midyind]]
            x, y  = radec2xy(RA, Dec, midRA, midDec)
        except IndexError:
            midRA = RA[0]
            midDec = Dec[0]
    else:
        midRA = RA[0]
        midDec = Dec[0]

    return np.array([x, y]), midRA, midDec


def radec2xy(RA, Dec, refRA=None, refDec=None):
    """
    Returns x, y for input ra, dec.

    Note that the reference RA and Dec must be the same in calls to both
    radec2xy() and xy2radec() if matched pairs of (x, y) <=> (RA, Dec) are
    desired.

    Parameters
    ----------
    RA : list
        List of RA values in degrees
    Dec : list
        List of Dec values in degrees
    refRA : float, optional
        Reference RA in degrees.
    refDec : float, optional
        Reference Dec in degrees

    Returns
    -------
    x, y : list, list
        Lists of x and y pixel values corresponding to the input RA and Dec
        values

    """
    import numpy as np

    x = []
    y = []
    if refRA is None:
        refRA = RA[0]
    if refDec is None:
        refDec = Dec[0]

    # Make wcs object to handle transformation from ra and dec to pixel coords.
    w = makeWCS(refRA, refDec)

    for ra_deg, dec_deg in zip(RA, Dec):
        ra_dec = np.array([[ra_deg, dec_deg]])
        x.append(w.wcs_world2pix(ra_dec, 0)[0][0])
        y.append(w.wcs_world2pix(ra_dec, 0)[0][1])

    return x, y


def xy2radec(x, y, refRA=0.0, refDec=0.0):
    """
    Returns x, y for input ra, dec.

    Note that the reference RA and Dec must be the same in calls to both
    radec2xy() and xy2radec() if matched pairs of (x, y) <=> (RA, Dec) are
    desired.

    Parameters
    ----------
    x : list
        List of x values in pixels
    y : list
        List of y values in pixels
    refRA : float, optional
        Reference RA in degrees
    refDec : float, optional
        Reference Dec in degrees

    Returns
    -------
    RA, Dec : list, list
        Lists of RA and Dec values corresponding to the input x and y pixel
        values

    """
    import numpy as np

    RA = []
    Dec = []

    # Make wcs object to handle transformation from ra and dec to pixel coords.
    w = makeWCS(refRA, refDec)

    for xp, yp in zip(x, y):
        x_y = np.array([[xp, yp]])
        RA.append(w.wcs_pix2world(x_y, 0)[0][0])
        Dec.append(w.wcs_pix2world(x_y, 0)[0][1])

    return RA, Dec


def makeWCS(refRA, refDec):
    """
    Makes simple WCS object.

    Parameters
    ----------
    refRA : float
        Reference RA in degrees
    refDec : float
        Reference Dec in degrees

    Returns
    -------
    w : astropy.wcs.WCS object
        A simple TAN-projection WCS object for specified reference position

    """
    from astropy.wcs import WCS
    import numpy as np

    w = WCS(naxis=2)
    w.wcs.crpix = [1000, 1000]
    w.wcs.cdelt = np.array([-0.066667, 0.066667])
    w.wcs.crval = [refRA, refDec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.set_pv([(2, 1, 45.0)])

    return w


def calculateSeparation(ra1, dec1, ra2, dec2):
    """
    Returns angular separation between two coordinates (all in degrees).

    Parameters
    ----------
    ra1 : float or numpy array
        RA of coordinate 1 in degrees
    dec1 : float or numpy array
        Dec of coordinate 1 in degrees
    ra2 : float
        RA of coordinate 2 in degrees
    dec2 : float
        Dec of coordinate 2 in degrees

    Returns
    -------
    separation : astropy Angle or numpy array
        Angular separation in degrees

    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    coord1 = SkyCoord(ra1, dec1, unit=(u.degree, u.degree), frame='fk5')
    coord2 = SkyCoord(ra2, dec2, unit=(u.degree, u.degree), frame='fk5')

    return coord1.separation(coord2)


# The following taken from
# http://code.activestate.com/recipes/578381-a-point-in-polygon-program-sw-sloan-algorithm/
def _det(xvert, yvert):
    """
    Compute twice the area of the triangle defined by points with using
    determinant formula.

    Parameters
    ----------
    xvert : array
        A vector of nodal x-coords (array-like).
    yvert : array
        A vector of nodal y-coords (array-like).

    Returns
    -------
    area : float
        Twice the area of the triangle defined by the points.
        _det is positive if points define polygon in anticlockwise order.
        _det is negative if points define polygon in clockwise order.
        _det is zero if at least two of the points are concident or if
            all points are collinear.

    """
    xvert = np.asfarray(xvert)
    yvert = np.asfarray(yvert)
    x_prev = np.concatenate(([xvert[-1]], xvert[:-1]))
    y_prev = np.concatenate(([yvert[-1]], yvert[:-1]))
    return np.sum(yvert * x_prev - xvert * y_prev, axis=0)


class Polygon:
    """
    Polygon object.

    Parameters
    ----------
    x : array
        A sequence of nodal x-coords.
    y : array
        A sequence of nodal y-coords.

    """

    def __init__(self, x, y):
        if len(x) != len(y):
            raise IndexError('x and y must be equally sized.')
        self.x = np.asfarray(x)
        self.y = np.asfarray(y)
        # Closes the polygon if were open
        x1, y1 = x[0], y[0]
        xn, yn = x[-1], y[-1]
        if x1 != xn or y1 != yn:
            self.x = np.concatenate((self.x, [x1]))
            self.y = np.concatenate((self.y, [y1]))
        # Anti-clockwise coordinates
        if _det(self.x, self.y) < 0:
            self.x = self.x[::-1]
            self.y = self.y[::-1]

    def is_inside(self, xpoint, ypoint, smalld=1e-12):
        '''Check if point is inside a general polygon.

        Input parameters:

        xpoint -- The x-coord of the point to be tested.
        ypoint -- The y-coords of the point to be tested.
        smalld -- A small float number.

        xpoint and ypoint could be scalars or array-like sequences.

        Output parameters:

        mindst -- The distance from the point to the nearest point of the
                  polygon.
                  If mindst < 0 then point is outside the polygon.
                  If mindst = 0 then point in on a side of the polygon.
                  If mindst > 0 then point is inside the polygon.

        Notes:

        An improved version of the algorithm of Nordbeck and Rydstedt.

        REF: SLOAN, S.W. (1985): A point-in-polygon program. Adv. Eng.
             Software, Vol 7, No. 1, pp 45-47.

        '''
        xpoint = np.asfarray(xpoint)
        ypoint = np.asfarray(ypoint)
        # Scalar to array
        if xpoint.shape is tuple():
            xpoint = np.array([xpoint], dtype=float)
            ypoint = np.array([ypoint], dtype=float)
            scalar = True
        else:
            scalar = False
        # Check consistency
        if xpoint.shape != ypoint.shape:
            raise IndexError('x and y has different shapes')
        # If snear = True: Dist to nearest side < nearest vertex
        # If snear = False: Dist to nearest vertex < nearest side
        snear = np.ma.masked_all(xpoint.shape, dtype=bool)
        # Initialize arrays
        mindst = np.ones_like(xpoint, dtype=float) * np.inf
        j = np.ma.masked_all(xpoint.shape, dtype=int)
        x = self.x
        y = self.y
        n = len(x) - 1  # Number of sides/vertices defining the polygon
        # Loop over each side defining polygon
        for i in range(n):
            d = np.ones_like(xpoint, dtype=float) * np.inf
            # Start of side has coords (x1, y1)
            # End of side has coords (x2, y2)
            # Point has coords (xpoint, ypoint)
            x1 = x[i]
            y1 = y[i]
            x21 = x[i + 1] - x1
            y21 = y[i + 1] - y1
            x1p = x1 - xpoint
            y1p = y1 - ypoint
            # Points on infinite line defined by
            #     x = x1 + t * (x1 - x2)
            #     y = y1 + t * (y1 - y2)
            # where
            #     t = 0    at (x1, y1)
            #     t = 1    at (x2, y2)
            # Find where normal passing through (xpoint, ypoint) intersects
            # infinite line
            t = -(x1p * x21 + y1p * y21) / (x21 ** 2 + y21 ** 2)
            tlt0 = t < 0
            tle1 = (0 <= t) & (t <= 1)
            # Normal intersects side
            d[tle1] = ((x1p[tle1] + t[tle1] * x21) ** 2 +
                       (y1p[tle1] + t[tle1] * y21) ** 2)
            # Normal does not intersects side
            # Point is closest to vertex (x1, y1)
            # Compute square of distance to this vertex
            d[tlt0] = x1p[tlt0] ** 2 + y1p[tlt0] ** 2
            # Store distances
            mask = d < mindst
            mindst[mask] = d[mask]
            j[mask] = i
            # Point is closer to (x1, y1) than any other vertex or side
            snear[mask & tlt0] = False
            # Point is closer to this side than to any other side or vertex
            snear[mask & tle1] = True
        if np.ma.count(snear) != snear.size:
            raise IndexError('Error computing distances')
        mindst **= 0.5
        # Point is closer to its nearest vertex than its nearest side, check if
        # nearest vertex is concave.
        # If the nearest vertex is concave then point is inside the polygon,
        # else the point is outside the polygon.
        jo = j.copy()
        jo[j == 0] -= 1
        area = _det([x[j + 1], x[j], x[jo - 1]], [y[j + 1], y[j], y[jo - 1]])
        mindst[~snear] = np.copysign(mindst, area)[~snear]
        # Point is closer to its nearest side than to its nearest vertex, check
        # if point is to left or right of this side.
        # If point is to left of side it is inside polygon, else point is
        # outside polygon.
        area = _det([x[j], x[j + 1], xpoint], [y[j], y[j + 1], ypoint])
        mindst[snear] = np.copysign(mindst, area)[snear]
        # Point is on side of polygon
        mindst[np.fabs(mindst) < smalld] = 0
        # If input values were scalar then the output should be too
        if scalar:
            mindst = float(mindst)
        return mindst


def read_vertices(filename):
    """
    Returns facet vertices
    """
    import pickle

    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def mask_vertices(mask_im, vertices_file):
    """
    Modify the input image to exclude regions outside of the polygon
    """
    import pyrap.images as pim

    ma = mask_im.coordinates()
    new_im = pim.image('',shape=mask_im.shape(), coordsys=ma)
    bool_mask = pim.image('',shape=mask_im.shape(), coordsys=ma)

    img_type = mask_im.imagetype()
    data = mask_im.getdata()
    bool_data = np.ones(data.shape)

    vertices = read_vertices(vertices_file)
    RAverts = vertices[0]
    Decverts = vertices[1]
    xvert = []
    yvert = []
    for RAvert, Decvert in zip(RAverts, Decverts):
        pixels = mask_im.topixel([1, 1, Decvert*np.pi/180.0,
            RAvert*np.pi/180.0])
        xvert.append(pixels[2]) # x -> Dec
        yvert.append(pixels[3]) # y -> RA
    poly = Polygon(xvert, yvert)

    # Find distance to nearest poly edge and unmask those that
    # are outside the facet (dist < 0)
    masked_ind = np.indices(data[0, 0].shape)
    dist = poly.is_inside(masked_ind[0], masked_ind[1])
    outside_ind = np.where(dist < 0.0)
    if len(outside_ind[0]) > 0:
        data[0, 0, masked_ind[0][outside_ind], [masked_ind[1][outside_ind]]] = 0
        bool_data[0, 0, masked_ind[0][outside_ind], [masked_ind[1][outside_ind]]] = 0

    new_im.putdata(data)
    bool_mask.putdata(bool_data)

    return new_im, bool_mask


def find_nearest(direction1, directions):
    """
    Finds nearest direction to input direction

    Parameters
    ----------
    direction1 : Direction object
        Target direction for which nearset direction is to be found
    directions : list
        List of directions to search. Should not include the target direction

    Returns
    -------
    direction : Direction object
        Nearest direction

    """
    sep = []
    for direction2 in directions:
        sep.append(calculateSeparation(direction1.ra, direction1.dec,
                            direction2.ra, direction2.dec).value)

    return directions[np.argmin(sep)]

