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
    """
    points = []
    inc = pi * (3 - sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = sqrt(1 - y*y)
        phi = k * inc
        points.append([cos(phi)*r, y, sin(phi)*r])
    return points


def generate_circle_points(x, y, z, a, b, c, u, v, w, n):
    """
    Returns list of 3d coordinates of points on a circle using the
    Rodrigues rotation formula
    """
    points = []
    offset = 2 * pi / float(n)
    u_2 = u**2
    v_2 = v**2
    w_2 = w**2
    dst = u_2 + v_2 + w_2
    sqrtdst = sqrt(dst)
    uxvywz =  - u*x - v*y - w*z
    b_v = b*v
    c_w = c*w
    a_u = a*u
    one = (a * (v_2 + w_2) - u*(b_v + c_w + uxvywz))
    two = (b * (u_2 + w_2) - v*(a_u + c_w + uxvywz))
    tre = (c * (u_2 + v_2) - w*(a_u + b_v + uxvywz))
    onep = sqrtdst * (-c*v + b*w - w*y + v*z)
    twop = sqrtdst * ( c*u - a*w + w*x - u*z)
    trep = sqrtdst * (-b*u + a*v - v*x + u*y)
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
    """
    return ((part1['x'] - part2['x'])**2 +
            (part1['y'] - part2['y'])**2 +
            (part1['z'] - part2['z'])**2)


def distance(part1, part2):
    """
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

