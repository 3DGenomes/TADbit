"""
30 Oct 2013


"""

import sys

from warnings  import catch_warnings, simplefilter
from itertools import combinations
from math      import pi, sqrt, cos, sin, acos
from copy      import deepcopy

import numpy as np
from numpy.random import shuffle as np_shuffle
from scipy.stats  import skew, kurtosis, norm as sc_norm
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

from pytadbit.eqv_rms_drms import rmsdRMSD_wrapper
from pytadbit.consistency import consistency_wrapper
from pytadbit.utils.extraviews import tadbit_savefig

def generate_sphere_points(n=100):
    """
    Returns list of 3d coordinates of points on a sphere using the
    Golden Section Spiral algorithm.

    :param n: number of points in the sphere

    :returns a sphere of radius 1, centered in the origin
    """
    points = []
    inc = pi * (3 - sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = sqrt(1 - y*y)
        phi = k * inc
        # points.append(dict((('x', cos(phi) * r),('y', y),('z', sin(phi) * r))))
        points.append((cos(phi) * r, y, sin(phi) * r))
    return points


def get_center_of_mass(x, y, z, zeros):
    """
    get the center of mass of a given object with list of x, y, z coordinates
    """
    xm = ym = zm = 0.
    size = len(x)
    subsize  = 0
    for i in range(size):
        if not zeros[i]:
            continue
        subsize += 1
        xm += x[i]
        ym += y[i]
        zm += z[i]
    xm /= subsize
    ym /= subsize
    zm /= subsize
    return xm, ym, zm


def mass_center(x, y, z, zeros):
    """
    Transforms coordinates according to the center of mass

    :param x: list of x coordinates
    :param y: list of y coordinates
    :param z: list of z coordinates
    """
    xm, ym, zm = get_center_of_mass(x, y, z, zeros)
    for i in range(len(x)):
        x[i] -= xm
        y[i] -= ym
        z[i] -= zm


# def generate_circle_points(x, y, z, a, b, c, u, v, w, n):
#     """
#     Returns list of 3d coordinates of points on a circle using the
#     Rodrigues rotation formula.
#
#     see *Murray, G. (2013). Rotation About an Arbitrary Axis in 3 Dimensions*
#     for details
#
#     :param x: x coordinate of a point somewhere on the circle
#     :param y: y coordinate of a point somewhere on the circle
#     :param z: z coordinate of a point somewhere on the circle
#     :param a: x coordinate of the center
#     :param b: y coordinate of the center
#     :param c: z coordinate of the center
#     :param u: 1st element of a vector in the same plane as the circle
#     :param v: 2nd element of a vector in the same plane as the circle
#     :param w: 3rd element of a vector in the same plane as the circle
#     :param n: number of points in the circle
#
#     TODO: try simplification for a=b=c=0 (and do the translation in the main
#           function)
#     """
#     points = []
#     offset = 2 * pi / float(n)
#     u_2 = u**2
#     v_2 = v**2
#     w_2 = w**2
#     dst = u_2 + v_2 + w_2
#     sqrtdst = sqrt(dst)
#     uxvywz =  - u*x - v*y - w*z
#     b_v = b*v
#     c_w = c*w
#     a_u = a*u
#     one = (a * (v_2 + w_2) - u*(b_v + c_w + uxvywz))
#     two = (b * (u_2 + w_2) - v*(a_u + c_w + uxvywz))
#     tre = (c * (u_2 + v_2) - w*(a_u + b_v + uxvywz))
#     onep = sqrtdst * (-c*v + b*w - w*y + v*z)
#     twop = sqrtdst * ( c*u - a*w + w*x - u*z)
#     trep = sqrtdst * (-b*u + a*v - v*x + u*y)
#     for k in range(int(n)):
#         ang = k * offset
#         cosang = cos(ang)
#         dcosang = cosang * dst
#         sinang = sin(ang)
#         points.append([(one * (1 - cosang) + x * dcosang + onep * sinang) / dst,
#                        (two * (1 - cosang) + y * dcosang + twop * sinang) / dst,
#                        (tre * (1 - cosang) + z * dcosang + trep * sinang) / dst]
#                       )
#     return points


def rotate_among_y_axis(x, y, z, angle):
    """
    Rotate and object with a list of x, y, z coordinates among its center of
    mass
    """
    xj = []
    yj = []
    zj = []
    for xi, yi, zi in zip(*(x, y, z)):
        #dist = square_distance((xi, yi, zi), center_of_mass)
        xj.append(xi*cos(angle) + zi*sin(angle))
        yj.append(yi)
        zj.append(xi*-sin(angle)+zi*cos(angle))
    return xj, yj, zj


def find_angle_rotation_improve_x(x, y, z, center_of_mass):
    """
    Finds the rotation angle needed to face the longest edge of the molecule
    """
    # find most distant point from center of mass:
    coords = list(zip(*(x, y, z)))
    xdst, ydst, zdst = max(coords, key=lambda i: square_distance(i, center_of_mass))
    dist = distance((xdst, ydst, zdst), center_of_mass)
    angle = acos((-xdst**2 - (dist + sqrt(dist**2 - xdst**2))) /
                 (2 * dist**2) + 1)
    return angle


def generate_circle_points(x, y, z, u, v, w, n):
    """
    Returns list of 3d coordinates of points on a circle using the
    Rodrigues rotation formula.

    see *Murray, G. (2013). Rotation About an Arbitrary Axis in 3 Dimensions*
    for details

    :param x: x coordinate of a point somewhere on the circle
    :param y: y coordinate of a point somewhere on the circle
    :param z: z coordinate of a point somewhere on the circle
    :param a: x coordinate of the center
    :param b: y coordinate of the center
    :param c: z coordinate of the center
    :param u: 1st element of a vector in the same plane as the circle
    :param v: 2nd element of a vector in the same plane as the circle
    :param w: 3rd element of a vector in the same plane as the circle
    :param n: number of points in the circle

    TODO: try simplification for a=b=c=0 (and do the translation in the main
          function)
    """
    points = []
    offset = 2 * pi / float(n)
    u_2 = u**2
    v_2 = v**2
    w_2 = w**2
    dst = u_2 + v_2 + w_2
    sqrtdst = sqrt(dst)
    uxvywz =  - u*x - v*y - w*z
    one = (-u * (uxvywz))
    two = (-v * (uxvywz))
    tre = (-w * (uxvywz))
    onep = sqrtdst * (- w*y + v*z)
    twop = sqrtdst * (+ w*x - u*z)
    trep = sqrtdst * (- v*x + u*y)
    for k in range(int(n)):
        ang = k * offset
        cosang = cos(ang)
        dcosang = cosang * dst
        sinang = sin(ang)
        points.append([(one * (1 - cosang) + x * dcosang + onep * sinang) / dst,
                       (two * (1 - cosang) + y * dcosang + twop * sinang) / dst,
                       (tre * (1 - cosang) + z * dcosang + trep * sinang) / dst]
                      )
    return points


def square_distance(part1, part2):
    """
    Calculates the square distance between two particles.

    :param part1: coordinate (dict format with x, y, z keys)
    :param part2: coordinate (dict format with x, y, z keys)

    :returns: square distance between two points in space
    """
    return ((part1[0] - part2[0])**2 +
            (part1[1] - part2[1])**2 +
            (part1[2] - part2[2])**2)

def fast_square_distance(x1, y1, z1, x2, y2, z2):
    """
    Calculates the square distance between two coordinates.

    :param part1: coordinate (dict format with x, y, z keys)
    :param part2: coordinate (dict format with x, y, z keys)

    :returns: square distance between two points in space
    """
    return ((x1 - x2)**2 +
            (y1 - y2)**2 +
            (z1 - z2)**2)


def distance(part1, part2):
    """
    Calculates the distance between two particles.

    :param part1: coordinate in list format (x, y, z)
    :param part2: coordinate in list format (x, y, z)

    :returns: distance between two points in space
    """
    return sqrt((part1[0] - part2[0])**2 +
                (part1[1] - part2[1])**2 +
                (part1[2] - part2[2])**2)


def angle_between_3_points(point1, point2, point3):
    """
    Calculates the angle between 3 particles

    Given three particles A, B and C, the angle g (angle ACB, shown below):

    ::


                          A
                         /|
                        /i|
                      c/  |
                      /   |
                     /    |
                    B )g  |b
                     \    |
                      \   |
                      a\  |
                        \h|
                         \|
                          C

    is given by the theorem of Al-Kashi:

    .. math::

      b^2 = a^2 + c^2 - 2ac\cos(g)

    :param point1: list of 3 coordinate for x, y and z
    :param point2: list of 3 coordinate for x, y and z
    :param point3: list of 3 coordinate for x, y and z

    :returns: angle in radians

    """
    a = distance(point2, point3)
    c = distance(point1, point2)
    b = distance(point1, point3)

    try:
        g = acos((a**2 - b**2 + c**2) / (2 * a * c))
    except ValueError:
        g = 0.
    return g


def calc_consistency(models, nloci, zeros, dcutoff=200):
    combines = list(combinations(models, 2))
    parts = [0 for _ in range(nloci)]
    for pm in consistency_wrapper([model['x'] for model in models],
                                  [model['y'] for model in models],
                                  [model['z'] for model in models],
                                  zeros,
                                  nloci, dcutoff, list(range(len(models))),
                                  len(models)):
        for i, p in enumerate(pm):
            parts[i] += p
    return [float(p)/len(combines) * 100 for p in parts]


def calc_eqv_rmsd(models, beg, end, zeros, dcutoff=200, one=False, what='score',
                  normed=True):
    """
    Calculates the RMSD, dRMSD, the number of equivalent positions and a score
    combining these three measures. The measure are done between a group of
    models in a one against all manner.

    :param beg: start particle number of the region to compare
    :param end: end particle number of the region to compare
    :param zeros: list of True/False representing particles to skip
    :param 200 dcutoff: distance in nanometer from which it is considered
       that two particles are separated.
    :param 0.75 fact: Factor for equivalent positions
    :param False one: if True assumes that only two models are passed, and
       returns the rmsd of their comparison
    :param 'score' what: values to return. Can be one of 'score', 'rmsd',
       'drmsd' or 'eqv'
    :param True normed: normalize result by maximum value (only applies to rmsd
       and drmsd)

    :returns: a score of each pairwise comparison according to:

       .. math::

         score_i = eqvs_i \\times \\frac{dRMSD_i / max(dRMSD)}
                                         {RMSD_i / max(RMSD)}

       where :math:`eqvs_i` is the number of equivalent position for the ith
       pairwise model comparison.

    """
    what = what.lower()
    if not what in ['score', 'rmsd', 'drmsd', 'eqv']:
        raise NotImplementedError("Only 'score', 'rmsd', 'drmsd' or 'eqv' " +
                                  "features are available\n")
    # remove particles with zeros from calculation
    x = []
    y = []
    z = []
    for m in range(len(models)):
        x.append([models[m]['x'][i] for i in range(beg, end) if zeros[i]])
        y.append([models[m]['y'][i] for i in range(beg, end) if zeros[i]])
        z.append([models[m]['z'][i] for i in range(beg, end) if zeros[i]])
    zeros = tuple([True for _ in range(len(x[0]))])
    scores = rmsdRMSD_wrapper(x, y, z, zeros, len(zeros),
                              dcutoff, list(range(len(models))), len(models),
                              int(one), what, int(normed))
    return scores


def dihedral(a, b, c, d, e):
    """
    Calculates dihedral angle between 4 points in 3D (array with x,y,z)
    """
    v1 = getNormedVector(b - a)
    v2 = getNormedVector(b - c)
    v4 = getNormedVector(d - c)
    v3 = getNormedVector(c - e)
    v1v2 = np.cross(v1, v2)
    v3v4 = np.cross(v3, v4)
    sign = 1 if np.linalg.det([v2, v1v2, v3v4]) < 0 else -1
    angle = getAngle(v1v2, v3v4)
    return sign * angle


def getNormedVector(dif):
    return (dif) / np.linalg.norm(dif)


def getAngle(v1v2, v2v3):
    return np.rad2deg(
        np.arccos(np.dot(
            v1v2   / np.linalg.norm(v1v2),
            v2v3.T / np.linalg.norm(v2v3)))
        )


def build_mesh(xis, yis, zis, nloci, nump, radius, superradius, include_edges):
    """
    Main function for the calculation of the accessibility of a model.
    """
    superradius = superradius or 1
    # number of dots in a circle is dependent the ones in a sphere
    numc = sqrt(nump) * sqrt(pi)
    right_angle = pi / 2 - pi / numc
    # keeps the remaining of integer conversion, to correct
    remaining = int(100*(numc - int(numc)) + 0.5)
    c_count = 0
    # number of circles per sphere needed to get previous equality are
    # dependent of:
    fact = float(nump)/numc/(2*radius)
    # starts big loop
    points    = [] # stores the particle coordinates and,
                   # if include_edges is True, the edge segments
    subpoints = [] # store the coordinates of each dot in the mesh
    supersubpoints = [] # store the coordinates of each dot in the mesh
    positions = {} # a dict to get dots belonging to a given point
    sphere    = generate_sphere_points(nump)
    i = 0
    for i in range(nloci - 1):
        modelx   = xis[i]
        modely   = yis[i]
        modelz   = zis[i]
        modelx1  = xis[i+1]
        modely1  = yis[i+1]
        modelz1  = zis[i+1]
        if i < nloci - 2:
            modelx2  = xis[i+2]
            modely2  = yis[i+2]
            modelz2  = zis[i+2]
        if i:
            modelx_1 = xis[i-1]
            modely_1 = yis[i-1]
            modelz_1 = zis[i-1]
        point = [modelx, modely, modelz]
        points.append(point)
        # get minimum length from next particle to display the sphere dot
        adj1 = distance(point, [modelx1, modely1, modelz1])

        # find a vector orthogonal to the axe between particle i and i+1
        difx = modelx - modelx1
        dify = modely - modely1
        difz = modelz - modelz1
        try:
            orthox = 1.
            orthoy = 1.
            orthoz = -(difx + dify) / difz
            #normer = sqrt(orthox**2 + orthoy**2 + orthoz**2) / radius
            normer = sqrt(2. + orthoz**2)# / radius
        except ZeroDivisionError:
            try:
                orthox = 1.
                orthoy = -(difx + difz) / dify
                orthoz = 1.
                #normer = sqrt(orthox**2 + orthoy**2 + orthoz**2) / radius
                normer = sqrt(2. + orthoz**2)# / radius
            except ZeroDivisionError:
                try:
                    orthox = 1.
                    orthoy = -(difx + difz) / dify
                    orthoz = 1.
                    #normer = sqrt(orthox**2 + orthoy**2 + orthoz**2) / radius
                    normer = sqrt(2. + orthoz**2)# / radius
                except ZeroDivisionError:
                    orthox = 1.
                    orthoy = 1.
                    orthoz = 1.
        orthox /= normer
        orthoy /= normer
        orthoz /= normer
        # define the number of circle to draw in this section
        between = int(fact * adj1 + 0.5)
        try:
            stepx = difx / between
            stepy = dify / between
            stepz = difz / between
        except ZeroDivisionError:
            stepx = stepy = stepz = 0

        hyp1 = sqrt(adj1**2 + radius**2)
        # this is an attempt of correction for the integrity of dots
        # uses intercept theorem
        hyp1 = (hyp1 - hyp1 / (2 * (1 + between)))**2

        # get minimum length from prev particle to display the sphere dot
        if i:
            adj2 = distance(point, [modelx_1, modely_1, modelz_1])
            hyp2 = sqrt(adj2**2 + radius**2)
            # this is an attempt of correction for the integrity of dots
            hyp2 = (hyp2 - hyp2 / (2 * (1 + between)))**2

        # set sphere around each particle
        for xxx, yyy, zzz in sphere:
            thing = [xxx * radius + modelx,
                     yyy * radius + modely,
                     zzz * radius + modelz]
            # same for super mesh
            superthing = [xxx * superradius + modelx,
                          yyy * superradius + modely,
                          zzz * superradius + modelz]
            # only place mesh outside torsion angle
            if fast_square_distance(modelx1, modely1, modelz1,
                                    thing[0], thing[1], thing[2]) > hyp1:
                if not i:
                    subpoints.append(thing)
                    supersubpoints.append(superthing)
                elif fast_square_distance(modelx_1, modely_1, modelz_1,
                                          thing[0], thing[1], thing[2]) > hyp2:
                    subpoints.append(thing)
                    supersubpoints.append(superthing)
                else:
                    continue
                positions.setdefault(i, []).append(len(subpoints)-1)

        def _add_circle(k, ptx, pty, ptz):
            for spoint in generate_circle_points(
                orthox, orthoy, orthoz, difx ,dify, difz,
                # correction for integer of numc
                numc + (1 if c_count%100 < remaining else 0)):
                dot = [spoint[0] * radius + ptx,
                       spoint[1] * radius + pty,
                       spoint[2] * radius + ptz]
                superdot = [spoint[0] * superradius + ptx,
                            spoint[1] * superradius + pty,
                            spoint[2] * superradius + ptz]
                # check that dot in circle is not too close from next edge
                if i < nloci - 2:
                    hyp = distance((modelx1, modely1, modelz1), dot)
                    ang = angle_between_3_points(dot,
                                                 (modelx1, modely1, modelz1),
                                                 (modelx2, modely2, modelz2))
                    if ang < right_angle:
                        if sin(ang) * hyp < radius:
                            continue
                # check that dot in circle is not too close from previous edge
                if i:
                    hyp = distance((modelx, modely, modelz), dot)
                    ang = angle_between_3_points(dot,
                                                 (modelx, modely, modelz),
                                                 (modelx_1, modely_1, modelz_1))
                    if ang < right_angle:
                        if sin(ang) * hyp < radius:
                            continue
                # print 'here'
                subpoints.append([dot[0], dot[1], dot[2]])
                supersubpoints.append([superdot[0], superdot[1], superdot[2]])
                positions.setdefault(i + float(k)/between, []).append(
                    len(subpoints) - 1)

        # define slices
        for k in range(between - 1, 0, -1):
            point = [modelx - k * stepx, modely - k * stepy, modelz - k * stepz]
            points.append(point)
            pointx, pointy, pointz = point

            if not include_edges:
                continue
            # define circles
            _add_circle(k, pointx, pointy, pointz)
            c_count += 1

    # add last point!!
    point = [xis[i+1], yis[i+1], zis[i+1]]
    points.append(point)
    # and its sphere
    adj = distance(point, [modelx, modely, modelz])
    hyp2 = sqrt(adj**2 + radius**2)
    hyp2 = (hyp2 - hyp2 / (2 * (1 + between)))**2
    for xxx, yyy, zzz in sphere:
        thing = [xxx * radius + modelx1,
                 yyy * radius + modely1,
                 zzz * radius + modelz1]
        superthing = [xxx * superradius + modelx1,
                      yyy * superradius + modely1,
                      zzz * superradius + modelz1]
        if fast_square_distance(modelx, modely, modelz,
                                thing[0], thing[1], thing[2]) > hyp2:
            subpoints.append(thing)
            supersubpoints.append(superthing)
        positions.setdefault(i+1, []).append(len(subpoints)-1)

    return points, subpoints, supersubpoints, positions


def randomize_matrix(data, savefig=None):
    size = len(data)
    rand_data = deepcopy(data)
    for d in range(size):
        diag = list(zip(*[list(range(d, size)), list(range(size - d))]))
        rdiag = diag[:]
        np_shuffle(rdiag)
        for v in range(len(diag)):
            val = data[diag[v][0]][diag[v][1]]
            a, b = rdiag[v][0], rdiag[v][1]
            rand_data[b][a] = rand_data[a][b] = val
    if savefig:
        plt.subplot(211)
        plt.imshow(np.log2(data), interpolation='none')
        plt.subplot(212)
        plt.imshow(np.log2(rand_data), interpolation='none')
        plt.savefig(savefig, format='pdf')
        plt.close('all')
    return rand_data


def mmp_score(matrix, nrand=10, verbose=False, savefig=None):
    """
    :param matrix: list of lists
    :param 10 nrand: number of randomizations
    :param None savefig: path where to save figure

    :returns: 1- MMP score which ranges from 0 (bad) to 1 (good), and 2- the
       expected correlation of the contact matrices of the modeled chromatin
       with the original Hi-C data (plus the 3- lower and 4- upper values
       expected in 95% of the cases)
    """
    data = np.array([np.array([v for v in l]) for l in matrix])

    if verbose:
        sys.stdout.write('  - getting EigenVectors\n')
    egval, _ = np.linalg.eigh(data)
    # sort eigenvalues/vectors
    idx = (-egval).argsort()
    egval = egval[idx]

    regvals = []

    if verbose:
        sys.stdout.write('  - randomization\n')
    for i in range(int(nrand)):
        if verbose:
            sys.stdout.write('\r    ' + str(i + 1) + ' / ' + str(nrand))
            sys.stdout.flush()
        regval, _ = np.linalg.eigh(randomize_matrix(data))
        regval = [abs(j) for j in regval]
        regval.sort(reverse=True)
        regvals.append( regval)
    if verbose:
        sys.stdout.write('\n')
    regvals = list(zip(*regvals))
    rvmean = []
    for rv in regvals:
        rvmean.append(np.mean(rv))
    total = sum(rvmean)/100
    rvmean = [i/total for i in rvmean]

    err = []
    for rv in regvals:
        rvstd = np.std(rv/total)
        err.append(2 * rvstd)

    zdata = sorted(np.log2([data[i][j] for i in range(len(data))
                            for j in range(i, len(data)) if data[i][j]]))
    skewness = skew(zdata)
    kurtness = kurtosis(zdata)

    if savefig:
        _ = plt.figure(figsize=(14, 8))
        gs = gridspec.GridSpec(7, 5, wspace=0.5, hspace=1.5)
        ax1 = plt.subplot(gs[:   , 0:3])
        ax2 = plt.subplot(gs[1:5 , 3: ])
        ax3 = plt.subplot(gs[5:7 , 3: ])
        with np.errstate(divide='ignore'):
            img = ax2.imshow(np.log2(data), interpolation='none')
        plt.colorbar(img, ax=ax2)

        if savefig:
            ax2.set_title('Original matrix', size=12)
            ax2.tick_params(axis='both', which='major', labelsize=10)
        ax2.set_xlabel('Bin')
        ax2.set_ylabel('Bin')

        normfit = sc_norm.pdf(zdata, np.mean(zdata), np.std(zdata))
        _ = ax3.plot(zdata, normfit, ':o', color='grey', alpha=.4,
                     markersize=.5)
        ax3.tick_params(axis='both', which='major', labelsize=10)

        ax3.hist(zdata, bins=20, density=True, alpha=0.7, color='r')
        ax3.set_xlabel('Z-score')
        ax3.set_ylabel('Frequency')
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        rcParams['axes.axisbelow']  = True
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        rcParams['axes.axisbelow']  = True
        rcParams['axes.grid']       = True
        rcParams['grid.color']      = 'w'
        rcParams['grid.linestyle']  = '-'
        rcParams['grid.linewidth']  = 2
        # rcParams['grid.alpha']      = .3
        ax1.minorticks_on()
        ax1.grid(ls='-', color='w', alpha=.3, lw=2, which='major')
        ax1.grid(ls='-', b=True, color='w', alpha=.3, lw=1, which='minor')
        ax1.spines['top'].set_color('none')
        ax1.spines['right'].set_color('none')
        ax1.spines['bottom'].set_color('none')
        ax1.spines['left'].set_color('none')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.set_xscale('log')
        try:
            ax1.set_axis_bgcolor((.9,.9,.9))
        except AttributeError:
            ax1.set_facecolor((.9,.9,.9))

        ax1.errorbar(list(range(1, 1 + len(rvmean))), rvmean, yerr=err, ecolor='red',
                     color='orange', lw=2,
                     label='%s randomizations' % (nrand))

    total = sum(abs(egval)) / 100
    egval = np.array(sorted([e/total for e in abs(egval)], reverse=True))

    for i in range(len(rvmean)):
        if rvmean[i] + err[i] > egval[i]:
            break
    signifidx = i

    size = len(data)

    sev = sum(egval[:signifidx]-rvmean[:signifidx])

    if savefig:
        ax1.plot(list(range(1, 1 + len(rvmean))), egval,
                 color='green', lw=2, label='Observed data')

        ax1.fill_between(list(range(1, 1 + len(rvmean))), rvmean, egval,
                         where=(np.array(rvmean) + np.array(err))<egval,
                         facecolor='green', interpolate=True, alpha=0.2)
        ax1.fill_between(list(range(1, 1 + len(rvmean))), rvmean, egval,
                         where=(np.array(rvmean) + np.array(err))>egval,
                         facecolor='red'  , interpolate=True, alpha=0.2)
        with catch_warnings():
            simplefilter("ignore")
            ax1.set_xlim((0,len(rvmean)))
            ax1.set_ylim((0, max(max(rvmean), max(egval))))
        ax1.legend(frameon=False, loc='upper right', prop={'size': 10})
        ax1.set_xlabel('Log indexes of Eigenvalues')
        ax1.set_ylabel('Eigenvalues (percentage of total)')
        #plt.subplots_adjust(right=0.6)

        #img = Image.open(opts.outdir + '/matrix_small.png')
        #fig.figimage(img, 640, -160)

    minv = float(min([i for d in data for i in d if i])) / 2
    if minv == 0.5:
        minv = 1./(len(data)**2)

    mmp = -0.0002 * size + 0.0335 * skewness - 0.0229 * kurtness + 0.0069 * sev + 0.8126

    if verbose:
        sys.stdout.write('\n')
        sys.stdout.write('\n                       Results\n')
        sys.stdout.write('                       -------\n\n')


    if verbose:
        sys.stdout.write('                  MMP score: %.4f\n\n' % mmp)

    ex_a1, ex_b1 = [0.6975926,  0.2548171]
    supa1, supb1 = [0.69300732000423904, 0.29858572176099613]
    lowa1, lowb1 = [0.70217788900976075, 0.211048473299004]

    scc     = (mmp - ex_b1 ) / ex_a1
    scc_up1 = (mmp - supb1 ) / supa1
    scc_lw1 = (mmp - lowb1 ) / lowa1

    if verbose:
        sys.stdout.write(('  predicted dSCC is %.3f (%.3f-%.3f '
                          '68%% confidence)\n') % (scc , scc_up1 , scc_lw1 ))

    supa75, supb75 = [0.69230778430383244, 0.30526310790548261]
    lowa75, lowb75 = [0.70287742471016734, 0.20437108715451746]

    scc_up75 = (mmp - supb75 ) / supa75
    scc_lw75 = (mmp - lowb75 ) / lowa75

    if verbose:
        sys.stdout.write(('                        (%.3f-%.3f '
                          '75%% confidence)\n') % (scc_up75 , scc_lw75 ))

    supa2, supb2 = [0.68855373600821357, 0.34109720480765293]
    lowa2, lowb2 = [0.70663147300578644, 0.16853699025234709]

    scc_up2 = (mmp - supb2 ) / supa2
    scc_lw2 = (mmp - lowb2 ) / lowa2
    if verbose:
        sys.stdout.write(('                        (%.3f-%.3f '
                          '95%% confidence)\n') % (scc_up2 , scc_lw2 ))

    if savefig:
        # write the log
        log = ''
        log +=  '    1- Matrix size (number of eigenvalues): %s\n' % (len(egval))
        log +=  "    2- Skewness of the distribution: %0.3f\n" % (skewness)
        log +=  "    3- Kurtosis of the distribution: %0.3f\n" % (kurtness)
        log +=  "    4- Sum of differences signif EV real-rand: %0.3f\n\n" % (sev)
        plt.figtext(0.62, 0.77, log, size='small')
        log =  "MMP score: %.3f\n" % (mmp)
        log +=  "Predicted dSCC: %.3f (%.3f-%.3f at 95%% conf)\n" % (scc, scc_up2, scc_lw2)
        plt.figtext(0.61, 0.87, log, size=12)
        tadbit_savefig(savefig)
        plt.close('all')

    return mmp, scc, scc_up2 , scc_lw2
