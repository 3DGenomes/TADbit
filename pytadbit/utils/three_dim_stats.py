"""
30 Oct 2013


"""

from pytadbit.eqv_rms_drms import rmsdRMSD_wrapper
from pytadbit.consistency import consistency_wrapper
from itertools import combinations
import numpy as np
from math import pi, sqrt, cos, sin


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


def dihedral(a,b,c,d):
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

