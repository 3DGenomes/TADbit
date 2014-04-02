"""
06 Aug 2013


"""

from bisect    import bisect_left
from itertools import combinations
from math      import log, log10, exp
from warnings  import warn
import numpy as np


def newton_raphson (guess, contour, sq_length, jmax=100, xacc=1e-12):
    """
    Newton-Raphson method as defined in:
    http://www.maths.tcd.ie/~ryan/TeachingArchive/161/teaching/newton-raphson.c.html
    used to search for the persistence length of a given model.
    
    :param guess: starting value
    :param contour: of the model
    :param sq_length: square of the distance between first and last particle
    :param 100 jmax: number of tries
    :param 1e-12 xacc: precision of the optimization

    :returns: persistence length
    """
    for _ in xrange(1, jmax):
        contour_guess = contour / guess
        expcx = exp(-contour_guess) - 1
        # function
        fx = 2 * pow(guess, 2) * (contour_guess + expcx) - sq_length
        # derivative
        df = 2 * contour * expcx + 4 * guess * (expcx + contour_guess)

        dx = -fx / df
        if abs(dx) < xacc:
            # print ("found root after %d attempts, at %f\n" % (j, x))
            return guess
        guess += dx
    raise Exception("ERROR: exceeded max tries no root\n")


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


def nozero_log(values):
    # Set the virtual minimum of the matrix to half the non-null real minimum
    minv = float(min(values)) / 2
    if minv > 1:
        warn('WARNING: probable problem with normalization, check.\n')
        minv /= 2  # TODO: something better
    logminv = log(minv)
    for i in xrange(values):
        try:
            values[i] = log(values[i])
        except ZeroDivisionError:
            values[i] = logminv
    

def zscore(values):
    """
    Calculates the log10, Z-score of a given list of values.
    
    .. note::
    
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
    # get the log trasnform values
    nozero_log(values)
    mean_v = np.mean(values.values())
    std_v  = np.std (values.values())
    # replace values by z-score
    for i in values:
        values[i] = (values[i] - mean_v) / std_v


def calinski_harabasz(scores, clusters):
    """
    Implementation of the CH score [CalinskiHarabasz1974]_, that has shown to be
    one the most accurate way to compare clustering methods
    [MilliganCooper1985]_ [Tibshirani2001]_.

    The CH score is:

    .. math::

        CH(k) = \\frac{B(k) / (k-1)}{W(k)/(n - k)}

    Where :math:`B(k)` and :math:`W(k)` are the between and within cluster sums
    of squares, with :math:`k` clusters, and :math:`n` the total number of
    elements (models in this case).
   
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

