"""
06 Aug 2013


"""

from bisect   import bisect_left
from itertools import combinations
from pytadbit.eqv_rms_drms import rmsdRMSD_wrapper
from pytadbit.consistency import consistency_wrapper
from math     import log10
import numpy as np


class Interpolate(object):
    """
    simple linear interpolation
    """
    def __init__(self, x_list, y_list):
        for i, (x, y) in enumerate(zip(x_list, x_list[1:])):
            if y - x < 0:
                raise ValueError("x must be in strictly ascending")
            if y - x == 0 and i >= len(x_list)-2:
                x_list = x_list[:-1]
                y_list = y_list[:-1]
        if any(y - x <= 0 for x, y in zip(x_list, x_list[1:])):
            raise ValueError("x must be in strictly ascending")
        x_list = self.x_list = map(float, x_list)
        y_list = self.y_list = map(float, y_list)
        intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
        self.slopes = [(y2 - y1)/(x2 - x1) for x1, x2, y1, y2 in intervals]
        
    def __call__(self, x):
        i = bisect_left(self.x_list, x) - 1
        return self.y_list[i] + self.slopes[i] * (x - self.x_list[i])


def zscore(values, size):
    """
    _______________________/___
                          /
                         /
                        /
                       /
                      /
                     /
                    /
                   /
                  /
                 /
                /
               /
              /                     score
          ___/_________________________________
            /

    """
    # do not take into account the diagonal
    nop = dict([(i + size * i,  None) for i in xrange(size)])
    minv = min([v for v in values if v]) / 2
    # get the log10 of values
    vals = [log10(v) if v > 0 and not v in nop else log10(minv) for v in values]
    mean_v = np.mean(vals)
    std_v = np.std(vals)
    # replace values by z-score
    for i in xrange(len(values)):
        if values[i] > 0:
            values[i] = (vals[i] - mean_v) / std_v
        elif values[i] == 0:
            values[i] = (log10(minv) - mean_v) / std_v
        else:
            values[i] = -99


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


def calinski_harabasz(scores, clusters):
    """
    Implementation of the CH score [CalinskiHarabasz1974]_, that has shown to be
    one the most accurate way to compare clustering methods
    [MilliganCooper1985]_ [Tibshirani2001]_.

    The CH score is:

    .. math::

        CH(k) = \\frac{B(k) / (k-1)}{W(k)/(n - k)}

    Where :math:`B(k)` and :math:`W(k)` are between and within cluster sums of
    squares, with :math:`k` clusters, and :math:`n` the total number of
    points (models in this case).
   
    :param scores: a dict with, as keys, a tuple with a pair of models; and, as
       value, the distance between these models.
    :param clusters: a dict with, as key, the cluster number, and as value a
       list of models
    :param nmodels: total number of models

    :returns: the CH score
    """
    cluster_list = [c for c in clusters if len(clusters[c]) > 1]
    if len(cluster_list) <= 1:
        return 0.0
    nmodels = sum([len(clusters[c]) for c in cluster_list])

    between_cluster = (sum([sum([sum([scores[(md1, md2)]**2
                                      for md1 in clusters[cl1]])
                                 for md2 in clusters[cl2]])
                            / (len(clusters[cl1]) * len(clusters[cl2]))
                            for cl1, cl2 in combinations(cluster_list, 2)])
                       / ((len(cluster_list) - 1.0) / 2))
    
    within_cluster = (sum([sum([scores[(md1, md2)]**2
                                for md1, md2 in combinations(clusters[cls], 2)])
                           / (len(clusters[cls]) * (len(clusters[cls]) - 1.0) / 2)
                           for cls in cluster_list]))
    
    return ((between_cluster / (len(cluster_list) - 1))
            /
            (within_cluster / (nmodels - len(cluster_list))))


def calc_eqv_rmsd(models, nloci, dcutoff=200, var='score', one=False):
    """
    :param nloci: number of particles per model
    :param 200 dcutoff: distance in nanometer from which it is considered
       that two particles are separated.
    :param 0.75 fact: Factor for equivalent positions
    :param 'score' var: value to return, can be either (i) 'drmsd' (symmetry
       independent: mirrors will show no differences) (ii) 'score' that is:

       ::

                               dRMSD[i] / max(dRMSD)
         score[i] = eqvs[i] * -----------------------
                                RMSD[i] / max(RMSD)

       where eqvs[i] is the number of equivalent position for the ith
       pairwise model comparison.
                                           
    :returns: a score (depends on 'var' argument)
    """
    scores = rmsdRMSD_wrapper([models[m]['x'] for m in xrange(len(models))],
                              [models[m]['y'] for m in xrange(len(models))],
                              [models[m]['z'] for m in xrange(len(models))],
                              nloci, dcutoff, range(len(models)),
                              len(models), int(one))
    # scores = {}
    # nrmsds = []
    # drmsds = []
    # for md1 in xrange(len(models)):
    #     md1s = models[md1]
    #     for md2 in xrange(md1 + 1, len(models)):
    #         md2s = models[md2]
    #         eqv, nrmsd, drmsd = rmsdRMSD_wrapper(
    #             md1s['x'], md1s['y'], md1s['z'],
    #             md2s['x'], md2s['y'], md2s['z'], nloci, dcutoff, 0)
    #         nrmsds.append(nrmsd)
    #         drmsds.append(drmsd)
    #         scores[(md1, md2)] = eqv * drmsd / nrmsd
    # if one:
    #     return drmsd
    # max_rmsd_ov_max_drmsd = max(nrmsds) / max(drmsds)
    # if var=='score':
    #     for md1, md2 in scores.keys()[:]:
    #         score = scores[(md1, md2)] * max_rmsd_ov_max_drmsd
    #         scores[(md1, md2)] = score
    #         scores[(md2, md1)] = score
    # elif var=='drmsd':
    #     for i, (md1, md2) in enumerate(scores.keys()):
    #         scores[(md2, md1)] = drmsds[i]
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


# from numpy import zeros
# from math import acos

# def crossproduct(ar1, ar2, base=1) :
#   """Calculate the cross-product of two vectors.

#   Positional arguments :
#   ar1    -- first vector (ndarray)
#   ar2    -- second vector (ndarray)
#           their base is given by the keyword argument of the same name

#   Keyword arguments :
#   base -- base index (default 1)

#   """
#   for ar_ in (ar1, ar2) :
#     if ar_ is None or 3 + base != len(ar_) :
#       raise Exception('Invalid parameter passed')

#   q1_ = ar1[base+1] * ar2[base+2] - ar1[base+2] * ar2[base+1]
#   q2_ = ar1[base+2] * ar2[base+0] - ar1[base+0] * ar2[base+2]
#   q3_ = ar1[base+0] * ar2[base+1] - ar1[base+1] * ar2[base+0]

#   ans = zeros(base + 3, 'd')
#   ans[base:] = (q1_, q2_, q3_)
  
#   return ans

# def spatproduct(ar1, ar2, ar3, base=1) :
#   """Calculate the scalar triple product of three vectors.

#   Positional arguments :
#   ar1    -- first vector (ndarray)
#   ar2    -- second vector (ndarray)
#   ar3    -- third vector (ndarray)
#           their base is given by the keyword argument of the same name

#   Keyword arguments :
#   base -- base index (default 1)

#   """
#   for ar_ in (ar1, ar2, ar3) :
#     if ar_ is None or 3 + base != len(ar_) :
#       raise Exception('Invalid parameter passed')

#   mat_ = zeros((3, 3), 'd')
#   mat_[0] = ar1[base:]
#   mat_[1] = ar2[base:]
#   mat_[2] = ar3[base:]

#   return np.linalg.det(mat_)

# def angle_vectors(ar1, ar2, base=1) :
#   """Calculate the angle between two vectors.

#   Positional arguments :
#   ar1    -- first vector (ndarray)
#   br2    -- second vector (ndarray)
#             their base is given by the keyword argument of the same name

#   Keyword arguments :
#   base -- base index (default 1)

#   Return the angle in grad.
  
#   """
#   for ar_ in (ar1, ar2) :
#     if ar_ is None or 3 + base != len(ar_) :
#       raise Exception('Invalid parameter passed')

#   ar1 = np.array(ar1[base:], 'd')
#   ar2 = np.array(ar2[base:], 'd')
    
#   scalar_prod = np.dot(ar1, ar2)
  
#   if 0. == scalar_prod :
#     return 90./np.pi

#   arg_ = scalar_prod / np.sqrt(np.dot(ar1, ar1) * np.dot(ar2, ar2))

#   # sometimes acos behaves strange...
#   if 1. < arg_ :
#     arg_ = 1.
#   elif -1. > arg_ :
#     arg_ = -1

#   return 180./np.pi * acos(arg_)

# def distance(ar1, ar2, base=1) :
#   """Calculate the distance between two points.

#   Positional arguments :
#   ar1    -- first point (ndarray)
#   ar2    -- second point (ndarray)
#             their base is given by the keyword argument of the same name

#   Keyword arguments :
#   base -- base index (default 1)
  
#   """
#   for ar_ in (ar1, ar2) :
#     if ar_ is None or 3 + base != len(ar_) :
#       raise Exception('Invalid parameter passed')
    
#   dist = np.array(ar1[base:], 'd') - np.array(ar2[base:], 'd')

#   return np.linalg.norm(dist)

# def angle(ar1, ar2, ar3, base=1) :
#   """Calculate the angle between three points.

#   Positional arguments :
#   ar1    -- first point (ndarray)
#   ar2    -- second point (ndarray)
#   ar3    -- third point (ndarray)
#             their base is given by the keyword argument of the same name

#   Keyword arguments :
#   base -- base index (default 1)

#   Return the angle in grad.
  
#   """
#   for ar_ in (ar1, ar2, ar3) :
#     if ar_ is None or 3 + base != len(ar_) :
#       raise Exception('Invalid parameter passed')
  
#   v_ba = np.array(ar1, 'd') - np.array(ar2, 'd')
#   v_bc = np.array(ar3, 'd') - np.array(ar2, 'd')

#   return angle_vectors(v_ba, v_bc, base)

# def dihedral_bis(ar1, ar2, ar3, ar4, base=0) :
#   """Calculate the dihedral angle between four points.

#   Positional arguments :
#   ar1    -- first point (ndarray)
#   ar2    -- second point (ndarray)
#   ar3    -- third point (ndarray)
#   ar4    -- fourth point (ndarray)
#             their base is given by the keyword argument of the same name

#   Keyword arguments :
#   base -- base index (default 1)

#   Return the dihedral angle in grad.
  
#   """
#   # for ar_ in (ar1, ar2, ar3, ar4) :
#   #   if ar_ is None or 3 + base != len(ar_) :
#   #     raise Exception('Invalid parameter passed')

#   v_ba = np.array(ar1, 'd') - np.array(ar2, 'd')
#   v_bc = np.array(ar3, 'd') - np.array(ar2, 'd')
#   v_cb = np.array(ar2, 'd') - np.array(ar3, 'd')
#   v_cd = np.array(ar4, 'd') - np.array(ar3, 'd')

#   # normal to the plane (a, b, c)
#   norm1 = crossproduct(v_ba, v_bc, base)

#   # normal to the plane (b, c, d)
#   norm2 = crossproduct(v_cb, v_cd, base)

#   # scalar triple product which defines the sign of the dihedral angle
#   if 0. > spatproduct(v_bc, norm1, norm2, base) :
#     sign = -1.
#   else :
#     sign = +1.
  
#   return sign * angle_vectors(norm1, norm2, base)
