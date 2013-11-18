"""
30 Oct 2013


"""

from pytadbit.eqv_rms_drms import rmsdRMSD_wrapper
from pytadbit.consistency import consistency_wrapper
from itertools import combinations
import numpy as np
from math import pi, sqrt, cos, sin, acos


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
    :param part1: coordinate (dict format with x, y, z keys)
    :param part2: coordinate (dict format with x, y, z keys)

    :returns: square distance between two points in space
    """
    return ((part1[0] - part2[0])**2 +
            (part1[1] - part2[1])**2 +
            (part1[2] - part2[2])**2)

def fast_square_distance(x1, y1, z1, x2, y2, z2):
    """
    :param part1: coordinate (dict format with x, y, z keys)
    :param part2: coordinate (dict format with x, y, z keys)

    :returns: square distance between two points in space
    """
    return ((x1 - x2)**2 +
            (y1 - y2)**2 +
            (z1 - z2)**2)


def distance(part1, part2):
    """
    :param part1: coordinate in list format (x, y, z)
    :param part2: coordinate in list format (x, y, z)

    :returns: distance between two points in space
    """
    return sqrt((part1[0] - part2[0])**2 +
                (part1[1] - part2[1])**2 +
                (part1[2] - part2[2])**2)


def angle_between_3_points(point1, point2, point3):
    """
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


def calc_consistency(models, nloci, dcutoff=200):
    combines = list(combinations(models, 2))
    parts = [0 for _ in xrange(nloci)]
    for pm in consistency_wrapper([models[m]['x'] for m in models],
                                  [models[m]['y'] for m in models],
                                  [models[m]['z'] for m in models],
                                  nloci, dcutoff, range(len(models)),
                                  len(models)):
        for i, p in enumerate(pm):
            parts[i] += p
    return [float(p)/len(combines) * 100 for p in parts]


def calc_eqv_rmsd(models, nloci, dcutoff=200, one=False):
    """

    :param nloci: number of particles per model
    :param 200 dcutoff: distance in nanometer from which it is considered
       that two particles are separated.
    :param 0.75 fact: Factor for equivalent positions
    :param False one: if True assumes that only two models are passed, and
       returns the rmsd of their comparison

    :returns: a score of each pairwise comparison according to:

       ::

                               dRMSD[i] / max(dRMSD)
         score[i] = eqvs[i] * -----------------------
                                RMSD[i] / max(RMSD)

       where eqvs[i] is the number of equivalent position for the ith
       pairwise model comparison.
       
    """
    scores = rmsdRMSD_wrapper([models[m]['x'] for m in xrange(len(models))],
                              [models[m]['y'] for m in xrange(len(models))],
                              [models[m]['z'] for m in xrange(len(models))],
                              nloci, dcutoff, range(len(models)),
                              len(models), int(one))
    return scores


def dihedral(a, b, c, d):
    """
    Calculates dihedral angle between 4 points in 3D (array with x,y,z)
    """
    v1 = getNormedVector(b - a)
    v2 = getNormedVector(b - c)
    v3 = getNormedVector(c - b)
    v4 = getNormedVector(c - d)
    v1v2 = np.cross(v1, v2)
    v2v3 = np.cross(v3, v4)
    sign = -1 if np.linalg.det([v2, v1v2, v2v3]) < 0 else 1
    angle = getAngle(v1v2, v2v3)
    angle, sign = (angle, sign) if angle <= 90 else (180 - angle, - sign)
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
    for i in xrange(nloci-1):
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
        point = dict((('x', modelx),
                      ('y', modely),
                      ('z', modelz)))
        point = [xis[i], yis[i], zis[i]]
        points.append(point)
        # get minimum length from next particle to display the sphere dot
        adj1 = distance(point, [modelx1, modely1, modelz1])

        # find a vector orthogonal to the axe between particle i and i+1
        difx = modelx - modelx1
        dify = modely - modely1
        difz = modelz - modelz1
        # orthox = 1.
        # orthoy = 1.
        orthoz = -(difx + dify) / difz
        #normer = sqrt(orthox**2 + orthoy**2 + orthoz**2) / radius
        normer = sqrt(2. + orthoz**2)# / radius
        orthox = 1. / normer
        orthoy = 1. / normer
        orthoz /= normer
        # define the number of circle to draw in this section
        between = int(fact * adj1 + 0.5)
        stepx = difx / between
        stepy = dify / between
        stepz = difz / between

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

        def _add_circle(k):
            for spoint in generate_circle_points(
                orthox, orthoy, orthoz, difx ,dify, difz,
                # correction for integer of numc
                numc + (1 if c_count%100 < remaining else 0)):
                dot = [spoint[0] * radius + pointx, spoint[1] * radius + pointy,
                       spoint[2] * radius + pointz]
                superdot = [spoint[0] * superradius + pointx, spoint[1] * superradius + pointy,
                            spoint[2] * superradius + pointz]
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
        for k in xrange(between - 1, 0, -1):
            point = [modelx - k * stepx, modely - k * stepy, modelz - k * stepz]
            points.append(point)
            pointx, pointy, pointz = point

            if not include_edges:
                continue
            # define circles
            _add_circle(k)
            c_count += 1

    # add last AND least point!!
    points.append([modelx1, modely1, modelz1])
    # and its sphere
    adj = distance(points[-1], [modelx, modely, modelz])
    hyp2 = sqrt(adj**2 + radius**2)
    hyp2 = (hyp2 - hyp2 * adj / (2 * between) / adj)**2
    for xxx, yyy, zzz in sphere:
        thing = [xxx * radius + modelx1,
                 yyy * radius + modely1,
                 zzz * radius + modelz1]
        superthing = [xxx * superradius + modelx,
                      yyy * superradius + modely,
                      zzz * superradius + modelz]
        if fast_square_distance(modelx1, modely1, modelz1,
                                thing[0], thing[1], thing[2]) > hyp2:
            subpoints.append(thing)
            supersubpoints.append(superthing)
        positions.setdefault(i+1, []).append(len(subpoints)-1)

    return points, subpoints, supersubpoints, positions



