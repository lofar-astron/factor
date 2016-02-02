"""
Module that holds all direction-related functions
"""
import os
import numpy as np
import logging
from factor.lib.direction import Direction
from factor.lib.polygon import Polygon
import sys
from scipy.spatial import Delaunay


log = logging.getLogger('factor:directions')


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
        types = np.dtype({
            'names': ['name', 'radec', 'atrous_do', 'mscale_field_do',
               'cal_imsize', 'solint_p', 'solint_a', 'dynamic_range',
                'region_selfcal', 'region_field', 'peel_skymodel', 'outlier_source',
                'cal_size_deg', 'cal_flux_mjy'],
            'formats':['S255', 'S255', 'S5', 'S5', int, int, int, 'S2',
                'S255', 'S255', 'S255', 'S5', float, float]})
        directions = np.genfromtxt(directions_file, comments='#', dtype=types)
    except ValueError:
        types = np.dtype({
            'names': ['name', 'radec', 'atrous_do', 'mscale_field_do',
              'cal_imsize', 'solint_p', 'solint_a', 'dynamic_range',
                'region_selfcal', 'region_field', 'peel_skymodel', 'outlier_source'],
            'formats':['S255', 'S255', 'S5', 'S5', int, int, int, 'S2',
                'S255', 'S255', 'S255', 'S5']})
        directions = np.genfromtxt(directions_file, comments='#', dtype=types)

    data = []
    for direction in directions:
        RAstr, Decstr = direction['radec'].split(',')
        ra = Angle(RAstr).to('deg').value
        dec = Angle(Decstr).to('deg').value

        # Check coordinates
        if np.isnan(ra) or ra < 0 or ra > 360:
            log.error('RA %f is wrong for direction: %s. Ignoring direction.'
            	% (direction['radec'], direction['name']))
            continue
        if np.isnan(dec) or dec < -90 or dec > 90:
            log.error('DEC %f is wrong for direction: %s. Ignoring direction.'
                % (direction['radec'], direction['name']))
            continue

        # Check atrous_do (wavelet) setting
        if direction['atrous_do'].lower() == 'empty':
            atrous_do = None
        elif direction['atrous_do'].lower() == 'true':
            atrous_do = True
        else:
            atrous_do = False

        # Check mscale_field_do (multi-scale) setting
        if direction['mscale_field_do'].lower() == 'empty':
            mscale_field_do = None
        elif direction['mscale_field_do'].lower() == 'true':
            mscale_field_do = True
        else:
            mscale_field_do = False

        # Check outlier_source (peeling) setting
        if direction['outlier_source'].lower() == 'empty':
            outlier_source = None
        elif direction['outlier_source'].lower() == 'true':
            outlier_source = True
        else:
            outlier_source = False

        # Set defaults
        if direction['solint_a'] < 0:
            direction['solint_a'] = 0 # 0 => set internally
        if direction['solint_p'] < 0:
            direction['solint_p'] = 0 # 0 => set internally
        if len(direction) > 13:
            if direction['cal_size_deg'] < 0.0 or np.isnan(direction['cal_size_deg']):
                cal_size_deg = None
            else:
                cal_size_deg = direction['cal_size_deg']
            if np.isnan(direction['cal_flux_mjy']):
                cal_flux_jy = None
            else:
                cal_flux_jy = direction['cal_flux_mjy'] / 1000.0
        else:
            cal_size_deg = None
            cal_flux_jy = None

        data.append(Direction(direction['name'], ra, dec, atrous_do,
        	mscale_field_do, direction['cal_imsize'], direction['solint_p'],
        	direction['solint_a'], direction['dynamic_range'],
        	direction['region_selfcal'], direction['region_field'],
        	direction['peel_skymodel'], outlier_source, factor_working_dir,
        	False, cal_size_deg, cal_flux_jy))

    return data


def make_directions_file_from_skymodel(s, flux_min_Jy, size_max_arcmin,
    directions_separation_max_arcmin, directions_max_num=None,
    interactive=False):
    """
    Selects appropriate calibrators from sky models and makes the directions file

    Parameters
    ----------
    s : LSMTool SkyModel object
        Skymodel made by grouping clean components of dir-independent model
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
    directions_file = 'factor_directions.txt'
    ds9_directions_file = 'factor_directions_ds9.reg'

    # Filter larger patches
    sizes = s.getPatchSizes(units='arcmin', weight=True)
    s.select(sizes < size_max_arcmin, aggregate=True, force=True)
    if len(s) == 0:
        log.critical("No sources found that meet the specified max size criteria.")
        sys.exit(1)
    log.info('Found {0} sources with sizes below {1} '
        'arcmin'.format(len(s.getPatchNames()), size_max_arcmin))

    # Filter fainter patches
    s.select('I > {0} Jy'.format(flux_min_Jy), aggregate='sum', force=True)
    if len(s) == 0:
        log.critical("No sources found that meet the specified min flux criteria.")
        sys.exit(1)
    log.info('Found {0} sources with fluxes above {1} Jy'.format(
        len(s.getPatchNames()), flux_min_Jy))

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

    # Trim directions list to get directions_max_num of directions
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


def make_initial_skymodel(band, max_radius_deg=None):
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
    max_radius_deg : float, optional
        Maximum radius in degrees from the phase center within which to include
        sources. If None, it is set to the FWHM (i.e., a diameter of 2 * FWHM)

    """
    import lsmtool

    # Set LSMTool logging level
    lsmtool._logging.setLevel('debug')

    # Load sky model
    s = lsmtool.load(band.skymodel_dirindep)

    # Group clean components by thresholding after convolving model with
    # 1-arcmin beam
    s.group('threshold', FWHM='60.0 arcsec', root='facet')
    s.remove('Patch = patch_*', force=True) # Remove sources that did not threshold
    if len(s) == 0:
        log.critical("No sources found through thresholding.")
        sys.exit(1)
    log.info('Found {0} sources through thresholding'.format(
        len(s.getPatchNames())))

    # Filter out sources that lie outside of maximum specific radius from phase
    # center
    if not hasattr(band, 'fwhm_deg'):
        band.set_image_sizes()
    if max_radius_deg is None:
        max_radius_deg = band.fwhm_deg # means a diameter of 2 * FWHM

    log.info('Removing sources beyond a radius of {0} degrees (corresponding to '
        'a diameter of {1} * FWHM of the primary beam at {2} MHz)...'.format(
        max_radius_deg, round(2.0*max_radius_deg/band.fwhm_deg, 1), band.freq/1e6))

    dist = s.getDistance(band.ra, band.dec, byPatch=True)
    s.remove(dist > max_radius_deg, aggregate=True)

    return s


def group_directions(directions, n_per_grouping={'1':0}, allow_reordering=True):
    """
    Sorts directions into groups that can be selfcaled simultaneously

    Directions are grouped by flux and then optionally reordered to maximize
    the miniumum separation between sources in a group

    Parameters
    ----------
    directions : list of Direction objects
        List of input directions to group
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
    direction_groups = []
    if n_per_grouping == {'1': 0}:
        for d in directions:
            direction_groups.append([d])
        log.debug('Processing each direction in series')
    else:
        def find_min_separation(group):
            """
            Finds the minimum separation in degrees between sources in a group
            """
            if len(group) == 1:
                return 0.0

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

        # Reorganize groups in each grouping level to maximize the separation
        # between directions. The separation is calculated as the weighted
        # total separation between members of the group
        if allow_reordering:
            log.info('Reordering directions to obtain max separation...')
            direction_groups_orig = direction_groups[:]
            remaining_directions = directions[:]
            fluxes = [d.apparent_flux_mjy for d in remaining_directions]
            if None in fluxes:
                use_fluxes = False
            else:
                use_fluxes = True
            if len(direction_groups) > 1:
                for i, group in enumerate(direction_groups_orig):
                    d0 = remaining_directions[0]
                    new_group = [d0]
                    remaining_directions.remove(d0)
                    ndir = len(group)
                    wsep_prev = [0] * len(remaining_directions)
                    if ndir > 1:
                        for j in range(1, ndir):
                            if use_fluxes:
                                weights = [d.apparent_flux_mjy for d in
                                    remaining_directions]
                            else:
                                weights = [len(remaining_directions)-k for k in
                                    range(len(remaining_directions))]
                            sep = [calculateSeparation(d0.ra, d0.dec,
                                d.ra, d.dec) for d in remaining_directions]
                            wsep_new = []
                            for s, w, wsep in zip(sep, weights, wsep_prev):
                                wsep_new.append(s.value*w + wsep)
                            d1 = remaining_directions[np.argmax(wsep_new)]
                            new_group.append(d1)
                            remaining_directions.remove(d1)
                            wsep_prev.pop(np.argmax(wsep_new))
                            wsep_new.pop(np.argmax(wsep_new))
                            wsep_prev = [p+n for p, n in zip(wsep_prev, wsep_new)]
                    direction_groups[i] = new_group

        log.debug('Processing directions in the following groups:')
        for i, group in enumerate(direction_groups):
            log.debug('Group {0}: {1}'.format(i+1, [d.name for d in group]))

    return direction_groups


def thiessen(directions_list, field_ra_deg, field_dec_deg, faceting_radius_deg,
    s=None, check_edges=False, target_ra=None, target_dec=None,
    target_radius_arcmin=None, beam_ratio=None):
    """
    Generates and add thiessen polygons or patches to input directions

    Parameters
    ----------
    directions_list : list of Direction objects
        List of input directions
    field_ra_deg : float
        RA in degrees of field center
    field_dec_deg : float
        Dec in degrees of field center
    faceting_radius_deg : float
        Maximum radius within which faceting will be done. Direction objects
        with postions outside this radius will get small rectangular patches
        instead of thiessen polygons
    s : LSMTool SkyModel object, optional
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
    beam_ratio : float, optional
        Ratio of semi-major (N-S) axis to semi-minor (E-W) axis for the primary
        beam

    """
    import lsmtool
    import shapely.geometry
    from shapely.ops import cascaded_union
    from itertools import combinations
    from astropy.coordinates import Angle

    # Select directions inside FOV (here defined as ellipse given by
    # faceting_radius_deg and the mean elevation)
    faceting_radius_pix = faceting_radius_deg / 0.066667 # radius in pixels
    field_x, field_y = radec2xy([field_ra_deg], [field_dec_deg],
        refRA=field_ra_deg, refDec=field_dec_deg)

    fx = []
    fy = []
    for th in range(0, 360, 1):
        fx.append(faceting_radius_pix * np.cos(th * np.pi / 180.0) + field_x[0])
        fy.append(faceting_radius_pix * beam_ratio * np.sin(th * np.pi / 180.0) + field_y[0])
    fov_poly_tuple = tuple([(xp, yp) for xp, yp in zip(fx, fy)])
    fov_poly = Polygon(fx, fy)

    points, _, _ = getxy(directions_list, field_ra_deg, field_dec_deg)
    for x, y, d in zip(points[0], points[1], directions_list):
        dist = fov_poly.is_inside(x, y)
        if dist < 0.0:
            # Source is outside of FOV, so use simple rectangular patches
            d.is_patch = True

    # Now do the faceting (excluding the patches)
    directions_list_thiessen = [d for d in directions_list if not d.is_patch]
    points, _, _ = getxy(directions_list_thiessen, field_ra_deg, field_dec_deg)
    points = points.T

    # Generate array of outer points used to constrain the facets
    nouter = 64
    means = np.ones((nouter, 2)) * points.mean(axis=0)
    offsets = []
    angles = [np.pi/(nouter/2.0)*i for i in range(0, nouter)]
    for ang in angles:
        offsets.append([np.cos(ang), np.sin(ang)])

    # Generate initial facets
    radius = 5.0 * faceting_radius_deg / 0.066667 # radius in pixels
    scale_offsets = radius * np.array(offsets)
    outer_box = means + scale_offsets
    points_all = np.vstack([points, outer_box])
    tri = Delaunay(points_all)
    circumcenters = np.array([_circumcenter(tri.points[t])
                              for t in tri.vertices])
    thiessen_polys = [_thiessen_poly(tri, circumcenters, n)
                      for n in range(len(points_all) - nouter)]

    # Check for vertices that are very close to each other, as this gives problems
    # to the edge adjustment below
    for thiessen_poly in thiessen_polys:
        dup_ind = 0
        for i, (v1, v2) in enumerate(zip(thiessen_poly[:-1], thiessen_poly[1:])):
            if (approx_equal(v1[0], v2[0], rel=1e-6) and
                approx_equal(v1[1], v2[1], rel=1e-6)):
                thiessen_poly.pop(dup_ind)
                dup_ind -= 1
            dup_ind += 1

    # Clip the facets at FOV
    for i, thiessen_poly in enumerate(thiessen_polys):
        polyv = np.vstack(thiessen_poly)
        poly_tuple = tuple([(xp, yp) for xp, yp in zip(polyv[:, 0], polyv[:, 1])])
        p1 = shapely.geometry.Polygon(poly_tuple)
        p2 = shapely.geometry.Polygon(fov_poly_tuple)
        if p1.intersects(p2):
            p1 = p1.intersection(p2)
            xyverts = [np.array([xp, yp]) for xp, yp in
                zip(p1.exterior.coords.xy[0].tolist(),
                p1.exterior.coords.xy[1].tolist())]
            thiessen_polys[i] = xyverts

    # Check for sources near / on facet edges and adjust regions accordingly
    if check_edges:
        log.info('Adjusting facets to avoid sources...')
        RA, Dec = s.getPatchPositions(asArray=True)
        sx, sy = radec2xy(RA, Dec, refRA=field_ra_deg, refDec=field_dec_deg)
        sizes = s.getPatchSizes(units='degree').tolist()

        if target_ra is not None and target_dec is not None and target_radius_arcmin is not None:
            log.info('Including target ({0}, {1}) in facet adjustment'.format(
                target_ra, target_dec))
            tra = Angle(target_ra).to('deg').value
            tdec = Angle(target_dec).to('deg').value
            tx, ty = radec2xy([tra], [tdec], refRA=field_ra_deg, refDec=field_dec_deg)
            sx.extend(tx)
            sy.extend(ty)
            sizes.append(target_radius_arcmin*2.0/1.2/60.0)

        # Set minimum size to 2*FWHM of resolution of high-res image
        fwhm = 2.0 * 25.0 / 3600.0 # degrees
        sizes = [max(size, fwhm) for size in sizes]

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

    # Add the final facet and patch info to the directions
    for d in directions_list:
        # Make calibrator patch
        sx, sy = radec2xy([d.ra], [d.dec], refRA=field_ra_deg, refDec=field_dec_deg)
        patch_width = d.cal_imsize * d.cellsize_selfcal_deg / 0.066667 # size of patch in pixels
        x0 = sx[0] - patch_width / 2.0
        y0 = sy[0] - patch_width / 2.0
        selfcal_poly = [np.array([x0, y0]),
                np.array([x0, y0+patch_width]),
                np.array([x0+patch_width, y0+patch_width]),
                np.array([x0+patch_width, y0])]

        if d.is_patch:
            # For sources beyond max radius, set facet poly to calibrator poly
            add_facet_info(d, selfcal_poly, selfcal_poly, field_ra_deg, field_dec_deg)
        else:
            poly = thiessen_polys[directions_list_thiessen.index(d)]
            add_facet_info(d, selfcal_poly, poly, field_ra_deg, field_dec_deg)


def add_facet_info(d, selfcal_poly, facet_poly, midRA, midDec):
    """
    Convert facet polygon from x, y to RA, Dec and find width of facet and
    facet center

    """
    poly_cal = np.vstack([selfcal_poly, selfcal_poly[0]])
    ra_cal, dec_cal = xy2radec(poly_cal[:, 0], poly_cal[:, 1], midRA, midDec)
    thiessen_poly_deg_cal = [np.array(ra_cal[0: -1]), np.array(dec_cal[0: -1])]

    poly = np.vstack([facet_poly, facet_poly[0]])
    ra, dec = xy2radec(poly[:, 0], poly[:, 1], midRA, midDec)
    thiessen_poly_deg = [np.array(ra[0: -1]), np.array(dec[0: -1])]

    # Find size and centers of facet regions in degrees
    xmin = np.min(poly[:, 0])
    xmax = np.max(poly[:, 0])
    xmid = xmin + (xmax - xmin) / 2.0
    ymin = np.min(poly[:, 1])
    ymax = np.max(poly[:, 1])
    ymid = ymin + (ymax - ymin) / 2.0

    ra1, dec1 = xy2radec([xmin], [ymin], midRA, midDec)
    ra2, dec2 = xy2radec([xmax], [ymax], midRA, midDec)
    ra3, dec3 = xy2radec([xmax], [ymin], midRA, midDec)
    ra_center, dec_center = xy2radec([xmid], [ymid], midRA, midDec)
    ra_width_deg = calculateSeparation(ra1, dec1, ra3, dec3)
    dec_width_deg = calculateSeparation(ra3, dec3, ra2, dec2)
    width_deg = max(ra_width_deg.value, dec_width_deg.value)

    d.vertices = thiessen_poly_deg
    d.vertices_cal = thiessen_poly_deg_cal
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
        lines.append('box({0}, {1}, {2}", {2}") # text={{{3}}}\n'.
            format(direction.ra, direction.dec, imsize, direction.name))
        lines.append('box({0}, {1}, {2}", {2}")\n'.
            format(direction.ra, direction.dec, imsize_unmasked))

    with open(outputfile, 'wb') as f:
        f.writelines(lines)


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


def getxy(directions_list, midRA=None, midDec=None):
    """
    Returns array of projected x and y values.

    Parameters
    ----------
    directions_list : list
        List of direction objects
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees

    Returns
    -------
    x, y : numpy array, numpy array, float, float
        arrays of x and y values

    """
    if len(directions_list) == 0:
        return np.array([0, 0]), 0, 0

    RA = []
    Dec = []
    for direction in directions_list:
        RA.append(direction.ra)
        Dec.append(direction.dec)

    if midRA is None or midDec is None:
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

    x, y  = radec2xy(RA, Dec, refRA=midRA, refDec=midDec)

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


def read_vertices(filename):
    """
    Returns facet vertices

    Parameters
    ----------
    filename : str
        Filename of pickle file that contains direction dictionary with
        vertices

    Returns
    -------
    vertices : list
        List of facet vertices

    """
    import pickle

    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def mask_vertices(mask_im, vertices_file):
    """
    Modify the input image to exclude regions outside of the polygon

    Parameters
    ----------
    mask_im : pyrap.images image() object
        Mask image to modify
    vertices_file: str
        Filename of pickle file that contains direction dictionary with
        vertices

    Returns
    -------
    new_im:  pyrap.images image() object
        Modified mask image
    bool_mask :  pyrap.images image() object
        Modified mask image

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
    sep : float
        Separation in degrees

    """
    sep = []
    for direction2 in directions:
        sep.append(calculateSeparation(direction1.ra, direction1.dec,
                            direction2.ra, direction2.dec).value)

    return directions[np.argmin(sep)], sep[np.argmin(sep)]


def _float_approx_equal(x, y, tol=1e-18, rel=1e-7):
    if tol is rel is None:
        raise TypeError('cannot specify both absolute and relative errors are None')
    tests = []
    if tol is not None: tests.append(tol)
    if rel is not None: tests.append(rel*abs(x))
    assert tests
    return abs(x - y) <= max(tests)


def approx_equal(x, y, *args, **kwargs):
    """approx_equal(float1, float2[, tol=1e-18, rel=1e-7]) -> True|False
    approx_equal(obj1, obj2[, *args, **kwargs]) -> True|False

    Return True if x and y are approximately equal, otherwise False.

    If x and y are floats, return True if y is within either absolute error
    tol or relative error rel of x. You can disable either the absolute or
    relative check by passing None as tol or rel (but not both).

    For any other objects, x and y are checked in that order for a method
    __approx_equal__, and the result of that is returned as a bool. Any
    optional arguments are passed to the __approx_equal__ method.

    __approx_equal__ can return NotImplemented to signal that it doesn't know
    how to perform that specific comparison, in which case the other object is
    checked instead. If neither object have the method, or both defer by
    returning NotImplemented, approx_equal falls back on the same numeric
    comparison used for floats.

    >>> almost_equal(1.2345678, 1.2345677)
    True
    >>> almost_equal(1.234, 1.235)
    False

    """
    if not (type(x) is type(y) is float):
        # Skip checking for __approx_equal__ in the common case of two floats.
        methodname = '__approx_equal__'
        # Allow the objects to specify what they consider "approximately equal",
        # giving precedence to x. If either object has the appropriate method, we
        # pass on any optional arguments untouched.
        for a,b in ((x, y), (y, x)):
            try:
                method = getattr(a, methodname)
            except AttributeError:
                continue
            else:
                result = method(b, *args, **kwargs)
                if result is NotImplemented:
                    continue
                return bool(result)
    # If we get here without returning, then neither x nor y knows how to do an
    # approximate equal comparison (or are both floats). Fall back to a numeric
    # comparison.
    return _float_approx_equal(x, y, *args, **kwargs)

