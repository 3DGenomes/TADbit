"""
16 Dec 2012


"""
from bisect import bisect_left
from numpy import mean, std
from math import log10

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


class matrix(object):
    """
    simple object to store hi-c matrices. Only keeps half matrix
    values MUST be integers
    matrices MUST be symetric
    """
    def __init__(self, nums, size):
        """
        
        """
        self.values = reduce_matrix(nums, size)
        self.nrows  = size


    def rebuild_matrix(self):
        values = [[[] for _ in xrange(self.nrows)] for _ in xrange(self.nrows)]
        for i in xrange(self.nrows):
            for j in xrange(self.nrows):
                values[i][j] = self.values[i*(j-i)]
        return values


def reduce_matrix(nums, size):
    return tuple([nums[i*j] for i in xrange(size) for j in xrange(i, size)])


def zscore(values):
    vals = []
    for v in values:
        if v > 0:
            vals.append(log10(v))
    mean_v = mean(vals)
    std_v = std(vals)
    for i in xrange(len(values)):
        if values[i] > 0:
            values[i] = (log10(values[i])-mean_v)/std_v
        else:
            values[i] = -99


def nicer(res):
    """
    writes resolution number for human beings.
    """
    if not res % 1000000000:
        return str(res)[:-9] + 'Gb'
    if not res % 1000000:
        return str(res)[:-6] + 'Mb'
    if not res % 1000:
        return str(res)[:-3] + 'Kb'
    return str(res) + 'b'


COLOR = {None: '\033[31m',
         0   : '\033[34m',
         1   : '\033[34m',
         2   : '\033[34m',
         3   : '\033[36m',
         4   : '\033[0m' ,
         5   : '\033[1m' ,
         6   : '\033[33m',
         7   : '\033[33m',
         8   : '\033[35m',
         9   : '\033[35m',
         10  : '\033[31m'}


def colorize(string, num):
    """
    Colorize with ANSII colors a string for printing in shell. this acording to
    a given number between 0 and 10

    :param string: the string to colorize
    :param num: a number between 0 and 10 (if None, number will be equal to 10)

    :returns: the string 'decorated' with ANSII color code
    """
    return '{}{}\033[m'.format(COLOR[num], string)
