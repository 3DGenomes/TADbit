"""
16 Dec 2012


"""
from bisect import bisect_left
from numpy import mean, std
from math import log10

try:
    from matplotlib import pyplot as plt
except ImportError:
    pass


def tad_breaker(tads, cluster, resolution, bins=20, show_plot=False,
                title=None, max_tad_size=3000000, test=False):
    """
    Find out if TAD boundaries are breaking clusters of colocalizing
    (epi)genetic elements.
    
    :param tads:
    :param cluster: cluster of elements
    :param title: for plotting
    :param resolution: in bases
    :param 20 bins: number of bins to use
    :param False show_plot:
    :para False test: for significant difference between extremes and center

    :returns: list of counts for each bin
    """
    def arrange(beg, end, num):
        step  = int((end - beg) / num)
        for k in xrange(int(beg), int(end), step):
            yield k
        yield int(end)
    brk_set = []
    for t in xrange(len(tads['start'])):
        #if not tads['score'][t] > 5:
        #    continue
        if tads['end'][t] - tads['start'][t] > max_tad_size:
            continue
        brk_set.append([b for b in arrange(tads['start'][t] * resolution,
                                           tads['end'][t] * resolution, bins)])
    counts = []
    for b in range(bins + 1):
        count = 0
        for start, end in cluster:
            for brk in brk_set:
                if start < brk[b] < end:
                    count += 1
                    break
        counts.append(1-float(count)/len(cluster))
    print counts
    if show_plot:
        plt.bar(range(21), counts)
        plt.title(title)
        plt.show()
    return counts


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


def zscore(values, size):
    nop = dict([(i + size * i,  None) for i in xrange(size)])
    vals = []
    for v in values:
        if v > 0 and not v in nop:
            vals.append(log10(v))
        elif v == 0:
            vals.append(0.0)
    mean_v = mean(vals)
    std_v = std(vals)
    for i in xrange(len(values)):
        # if i in nop:
        #     values[i] = float('nan')
        #     continue
        if values[i] > 0:
            values[i] = (log10(values[i])-mean_v)/std_v
        elif values[i] == 0:
            values[i] = -1
        else:
            print 'hola', values[i]
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


COLOR = {None: '\033[31m', # red
         0   : '\033[34m', # blue
         1   : '\033[34m', # blue
         2   : '\033[34m', # blue
         3   : '\033[36m', # cyan
         4   : '\033[0m' , # white
         5   : '\033[1m' , # bold white
         6   : '\033[33m', # yellow
         7   : '\033[33m', # yellow
         8   : '\033[35m', # purple
         9   : '\033[35m', # purple
         10  : '\033[31m'  # red
         }

COLORHTML = {None: '<span style="color:red;">'       , # red
             0   : '<span>'                          , # blue
             1   : '<span style="color:blue;">'      , # blue
             2   : '<span style="color:blue;">'      , # blue
             3   : '<span style="color:purple;">'    , # purple
             4   : '<span style="color:purple;">'    , # purple
             5   : '<span style="color:teal;">'      , # cyan
             6   : '<span style="color:teal;">'      , # cyan
             7   : '<span style="color:olive;">'     , # yellow
             8   : '<span style="color:olive;">'     , # yellow
             9   : '<span style="color:red;">'       , # red
             10  : '<span style="color:red;">'         # red
             }


def colorize(string, num, ftype='ansi'):
    """
    Colorize with ANSII colors a string for printing in shell. this acording to
    a given number between 0 and 10

    :param string: the string to colorize
    :param num: a number between 0 and 10 (if None, number will be equal to 10)

    :returns: the string 'decorated' with ANSII color code
    """
    color = COLOR if ftype=='ansi' else COLORHTML
    return '{}{}{}'.format(color[num], string,
                           '\033[m' if ftype=='ansi' else '</span>')
