"""
16 Dec 2012


"""
from bisect   import bisect_left
from numpy    import std
from math     import log10
from re       import sub
from warnings import warn
import numpy as np

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')

###TEST to be removed


# DAVIDE!!!
# xnam = 'TR2'
# for crm in ['2L', '2R', '3L', '3R', '4', 'X']:
#     crmbit = load_chromosome('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/chr{0}.tdb'.format(crm))
#     exp = crmbit.experiments[xnam]
#     exp.load_experiment('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
#     matrix = exp.get_hic_matrix()
#     means = [np.mean(c) for c in matrix]
#     mean = np.mean(means)
#     stde = np.std(means)
#     plt.clf()
#     plt.figure(figsize=(12, 12))
#     hist = plt.hist(means, bins=100, alpha=.3, color='grey')
#     a = plt.vlines(mean-stde*1.25, 0, 110, colors='g', linestyles='dashed')
#     b = plt.vlines(mean+stde*1.25, 0, 110, colors='r', linestyles='dashed')
#     plt.legend([a, b], ['lower cut: mean of means minus (stdev of means) times 1.25',
#                         'upper cut: mean of means plus (stdev of means) times 1.25'])
#     plt.title(crm)
#     plt.xlabel('mean value of the rows/columns')
#     # plt.show()
#     plt.savefig('/home/fransua/Pictures/dmel_dbau_' + crm + '.pdf', filtype='pdf')
#
#
# xnam = 'k562'
# for crm in [str(c) for c in range(1,23)] + ['X']:
#     crmbit = load_chromosome('/home/fransua/db/hi-c/dekker_hsap/hsap/100Kb/{0}/chr{0}.tdb'.format(crm))
#     exp = crmbit.experiments[xnam]
#     exp.load_experiment('/home/fransua/db/hi-c/dekker_hsap/hsap/100Kb/{0}/{1}.matrix'.format(crm, xnam))
#
# xnam = 'mESC_rep1'
# for crm in [str(c) for c in range(1,20)] + ['X']:
#     crmbit = load_chromosome('/home/fransua/db/hi-c/dixon_hsap-mmus/mmus/100Kb/{0}/chr{0}.tdb'.format(crm))
#     exp = crmbit.experiments[xnam]
#     exp.load_experiment('/home/fransua/db/hi-c/dixon_hsap-mmus/mmus/20Kb/{0}/{1}.matrix'.format(crm, xnam), resolution=100000)
#     # crmbit.visualize('BR', paint_tads=True, show=True)
#     matrix = exp.get_hic_matrix()
#     cols = []
#     for c in sorted(matrix, key=lambda x: sum(x)):
#         cols.append(sum(c))
#        
#     plt.clf()
#     plt.figure(figsize=(12, 12))
#     median = np.median(cols)
#     mad = np.median([abs(median - c ) for c in cols])
#     best =(None, None, None)
#     hist = plt.hist(cols, bins=100, alpha=.3, color='grey')
#     xp = range(0, cols[-1])
#     x, y, _ = hist
#     while list(x).count(0) > len(x)/2:
#         cols = cols[:-1]
#         plt.clf()
#         hist = plt.hist(cols, bins=100, alpha=.3, color='grey')
#         xp = range(0, cols[-1])
#         x, y, _ = hist
#     for order in range(4, 18):
#         z = np.polyfit(y[1:], x, order)
#         # print np.roots(np.polyder(z))
#         if np.roots(np.polyder(z))[-1] <= 0:
#             continue
#         if np.roots(np.polyder(z))[-1] >= median:
#             continue
#         p  = np.poly1d(z)
#         R2 = get_r2(p, x, y[1:])
#         if best[0] < R2:
#             best = (R2, order, p, z)
#         print order, R2
#     print best[1]
#     p, z = best[2:]
#
#     a = plt.plot(xp, p(xp), "--", color='k')
#     plt.title(crm)
#     b = plt.vlines(np.roots(np.polyder(z))[-1], 0, 110, colors='r', linestyles='dashed')
#     c = plt.vlines(median - mad * 1.5, 0, 110, colors='g', linestyles='dashed')
#     plt.legend(a+[b]+[c], ['polyfit \n{}'.format(
#         ''.join([sub('e([-+][0-9]+)', 'e^{\\1}', '${}{:.1}x^{}$'.format('+' if j>0 else '',
#                                                                         j,
#                                                                         '{' + str(i) + '}'))
#                  for i, j in enumerate(list(p)[::-1])])),
#                            'first solution of polynomial derivation',
#                            'median - (1.5 * median absolute deviation)'],
#                fontsize='x-small')
#     # plt.bar(cols, range(2303))
#     # plt.show()
#    
#     plt.savefig('/home/fransua/Pictures/dixon_' + crm + '.pdf', filtype='pdf')

######

def get_r2 (fun, X, Y, *args):
    sstot = sum([(Y[i]-np.mean(Y))**2 for i in xrange(len(Y))])
    sserr = sum([(Y[i] - fun(X[i], *args))**2 for i in xrange(len(Y))])
    return 1 - sserr/sstot


def filter_with_polynomial_fit(matrx, draw_hist=False):
    """
    fits the distribution of Hi-C interaction count by column in the matrix to
    a polynomial. Then searches for the first possible 
    """
    nbins = 100
    # get sum of columns
    cols = []
    for c in sorted(matrx, key=sum):
        cols.append(sum(c))
    cols = np.array(cols)
    if draw_hist:
        plt.figure(figsize=(9, 9))
    median = np.median(cols)
    # mad = np.median([abs(median - c ) for c in cols])
    best =(None, None, None, None)
    # bin the sum of columns
    xmin = min(cols)
    xmax = max(cols)
    y = np.linspace(xmin, xmax, nbins)
    hist = np.digitize(cols, y)
    x = [sum(hist == i) for i in range(1, nbins + 1)]
    if draw_hist:
        hist = plt.hist(cols, bins=100, alpha=.3, color='grey')
    xp = range(0, cols[-1])
    # check if the binning is correct
    # we want at list half of the bins with some data
    while list(x).count(0) > len(x)/2:
        cols = cols[:-1]
        xmin = min(cols)
        xmax = max(cols)
        y = np.linspace(xmin, xmax, nbins)
        hist = np.digitize(cols, y)
        x = [sum(hist == i) for i in range(1, nbins + 1)]
        if draw_hist:
            plt.clf()
            hist = plt.hist(cols, bins=100, alpha=.3, color='grey')
        xp = range(0, cols[-1])
    # find best polynomial fit in a given range
    for order in range(4, 14):
        z = np.polyfit(y, x, order)
        zpp = np.roots(np.polyder(z, m=2))
        roots = np.roots(np.polyder(z))
        # check that we are concave down, otherwise take next root
        pente = np.polyval(zpp, abs(roots[-2] - roots[-1]) / 2 + roots[-1])
        if pente > 0:
            root = roots[-1]
        else:
            root = roots[-2]
        # root must be higher than zero
        if root <= 0:
            continue
        # and lower than the median
        if root >= median:
            continue
        p  = np.poly1d(z)
        R2 = get_r2(p, x, y)
        if best[0] < R2:
            best = (R2, order, p, z, root)
    p, z, root = best[2:]
    if draw_hist:
        a = plt.plot(xp, p(xp), "--", color='k')
        b = plt.vlines(np.roots(np.polyder(z))[-1], 0, 110,
                       colors='r', linestyles='dashed')
        # c = plt.vlines(median - mad * 1.5, 0, 110, colors='g',
        #                linestyles='dashed')
        plt.legend(a+[b], ['polyfit \n{}'.format(
            ''.join([sub('e([-+][0-9]+)', 'e^{\\1}',
                         '${}{:.1}x^{}$'.format('+' if j>0 else '', j,
                                                '{' + str(i) + '}'))
                     for i, j in enumerate(list(p)[::-1])])),
                               'first solution of polynomial derivation'],
                   fontsize='x-small')
        # plt.legend(a+[b]+[c], ['polyfit \n{}'.format(
        #     ''.join([sub('e([-+][0-9]+)', 'e^{\\1}',
        #                  '${}{:.1}x^{}$'.format('+' if j>0 else '', j,
        #                                         '{' + str(i) + '}'))
        #              for i, j in enumerate(list(p)[::-1])])),
        #                        'first solution of polynomial derivation',
        #                        'median - (1.5 * median absolute deviation)'],
        #            fontsize='x-small')
        plt.show()
    # label as bad the columns with sums lower than the root
    bads = []
    for i, col in enumerate(matrx):
        if sum(col) < root:
            bads.append(i)
    # now stored in Experiment._zeros, used for getting more accurate z-scores
    return bads


def filter_by_stdev(matrx):
    means = [np.mean(c) for c in matrx]
    mean = np.mean(means)
    stde = np.std(means)
    root = mean - stde * 1.25
    # label as bad the columns with sums lower than the root
    bads = []
    for i, col in enumerate(matrx):
        if sum(col) < root:
            bads.append(i)
    # now stored in Experiment._zeros, used for getting more accurate z-scores
    return bads


def filter_by_mad(matrx):
    # get sum of columns
    cols = []
    for c in sorted(matrx, key=sum):
        cols.append(sum(c))
    cols = np.array(cols)
    median = np.median(cols)
    mad = np.median([abs(median - c ) for c in cols])
    root = median - mad * 1.5
    # label as bad the columns with sums lower than the root
    bads = []
    for i, col in enumerate(matrx):
        if sum(col) < root:
            bads.append(i)
    # now stored in Experiment._zeros, used for getting more accurate z-scores
    return bads


def hic_filtering_for_modelling(matrx, method='polynomial'):
    """
    :param matrx: Hi-C matrix of a given experiment
    :param polynomial method: method to use for filtering Hi-C columns. Aims to
       remove columns with abnormally low count of interactions

    :returns: the indexes of the columns not to be considered for the
       calculation of the z-score
    """
    if method == 'polynomial':
        bads = filter_with_polynomial_fit(matrx)
    elif method == 'mad':
        bads = filter_by_mad(matrx)
    elif method == 'stdev':
        bads = filter_by_stdev(matrx)
    else:
        raise Exception
    return bads


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
    # do not take into acount the diagonal
    nop = dict([(i + size * i,  None) for i in xrange(size)])
    vals = []
    for v in values:
        if v > 0 and not v in nop:
            vals.append(log10(v))
        elif v == 0:
            vals.append(0.0)
    mean_v = np.mean(vals)
    std_v = std(vals)
    for i in xrange(len(values)):
        # if i in nop:
        #     values[i] = float('nan')
        #     continue
        if values[i] > 0:
            values[i] = (log10(values[i]) - mean_v) / std_v
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
