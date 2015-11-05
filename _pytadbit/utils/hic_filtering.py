"""
06 Aug 2013
"""

from warnings                  import warn
from sys                       import stderr
from re                        import sub
from pytadbit.utils.extraviews import tadbit_savefig

import numpy as np

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')


def get_r2 (fun, X, Y, *args):
    sstot = sum([(Y[i]-np.mean(Y))**2 for i in xrange(len(Y))])
    sserr = sum([(Y[i] - fun(X[i], *args))**2 for i in xrange(len(Y))])
    return 1 - sserr/sstot


def filter_by_mean(matrx, draw_hist=False, silent=False, bads=None, savefig=None):
    """
    fits the distribution of Hi-C interaction count by column in the matrix to
    a polynomial. Then searches for the first possible 
    """
    nbins = 100
    if not bads:
        bads = {}
    # get sum of columns
    cols = []
    size = len(matrx)
    for c in sorted([[matrx.get(i+j*size, 0) for j in xrange(size) if not j in bads]
                     for i in xrange(size) if not i in bads], key=sum):
        cols.append(sum(c))
    cols = np.array(cols)
    if draw_hist:
        plt.figure(figsize=(9, 9))
    try:
        percentile = np.percentile(cols, 5)
    except IndexError:
        warn('WARNING: no columns to filter out')
        return bads
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
    xp = range(0, int(cols[-1]))
    # check if the binning is correct
    # we want at list half of the bins with some data
    try:
        cnt = 0
        while list(x).count(0) > len(x)/2:
            cnt += 1
            cols = cols[:-1]
            xmin = min(cols)
            xmax = max(cols)
            y = np.linspace(xmin, xmax, nbins)
            hist = np.digitize(cols, y)
            x = [sum(hist == i) for i in range(1, nbins + 1)]
            if draw_hist:
                plt.clf()
                hist = plt.hist(cols, bins=100, alpha=.3, color='grey')
            xp = range(0, int(cols[-1]))
            if cnt > 10000:
                raise ValueError
        # find best polynomial fit in a given range
        for order in range(6, 18):
            z = np.polyfit(y, x, order)
            zp = np.polyder(z, m=1)
            roots = np.roots(np.polyder(z))
            # check that we are concave down, otherwise take next root
            pente = np.polyval(zp, abs(roots[-2] - roots[-1]) / 2 + roots[-1])
            if pente > 0:
                root = roots[-1]
            else:
                root = roots[-2]
            # root must be higher than zero
            if root <= 0:
                continue
            # and lower than the median
            if root >= percentile:
                continue
            p  = np.poly1d(z)
            R2 = get_r2(p, x, y)
            # try to avoid very large orders by weigthing negatively their fit
            if order > 13:
                R2 -= float(order)/30
            if best[0] < R2:
                best = (R2, order, p, z, root)
        try:
            p, z, root = best[2:]
            if draw_hist:
                xlims = plt.xlim()
                ylims = plt.ylim()
                a = plt.plot(xp, p(xp), "--", color='k')
                b = plt.vlines(root, 0, plt.ylim()[1], colors='r', linestyles='dashed')
                # c = plt.vlines(median - mad * 1.5, 0, 110, colors='g',
                #                linestyles='dashed')
                try:
                    plt.legend(a+[b], ['polyfit \n%s' % (
                        ''.join([sub('e([-+][0-9]+)', 'e^{\\1}',
                                     '$%s%.1fx^%s$' % ('+' if j>0 else '', j,
                                                       '{' + str(i) + '}'))
                                 for i, j in enumerate(list(p)[::-1])])),
                                           'first solution of polynomial derivation'],
                               fontsize='x-small')
                except TypeError:
                    plt.legend(a+[b], ['polyfit \n%s' % (
                        ''.join([sub('e([-+][0-9]+)', 'e^{\\1}',
                                     '$%s%.1fx^%s$' % ('+' if j>0 else '', j,
                                                       '{' + str(i) + '}'))
                                 for i, j in enumerate(list(p)[::-1])])),
                                           'first solution of polynomial derivation'])
                # plt.legend(a+[b]+[c], ['polyfit \n{}'.format (
                #     ''.join([sub('e([-+][0-9]+)', 'e^{\\1}',
                #                  '${}{:.1}x^{}$'.format ('+' if j>0 else '', j,
                #                                         '{' + str(i) + '}'))
                #              for i, j in enumerate(list(p)[::-1])])),
                #                        'first solution of polynomial derivation',
                #                        'median - (1.5 * median absolute deviation)'],
                #            fontsize='x-small')
                plt.ylim([0, ylims[1]])
                plt.xlim(xlims)
                plt.xlabel('Sum of interactions')
                plt.xlabel('Number of columns with a given value')
                if savefig:
                    tadbit_savefig(savefig)
                else:
                    plt.show()
            # label as bad the columns with sums lower than the root
            for i, col in enumerate([[matrx.get(i+j*size, 0)
                                      for j in xrange(size)]
                                     for i in xrange(size)]):
                if sum(col) < root:
                    bads[i] = sum(col)
            # now stored in Experiment._zeros, used for getting more accurate z-scores
            if bads and not silent:
                stderr.write(('\nWARNING: removing columns having less than %s ' +
                              'counts:\n %s\n') % (
                                 round(root, 3), ' '.join(
                                     ['%5s'%str(i + 1) + (''if (j + 1) % 20 else '\n')
                                      for j, i in enumerate(sorted(bads.keys()))])))
        except:
            if not silent:
                stderr.write('WARNING: Too many zeroes to filter columns.' +
                             ' SKIPPING...\n')
            if draw_hist:
                plt.xlabel('Sum of interactions')
                plt.xlabel('Number of columns with a given value')
                if savefig:
                    tadbit_savefig(savefig)
                else:
                    plt.show()
    except ValueError:
        if not silent:
            stderr.write('WARNING: Too few data to filter columns based on ' +
                         'mean value.\n')
    if draw_hist:
        plt.close('all')
    return bads


def filter_by_zero_count(matrx, perc_zero, silent=True):
    """
    :param matrx: Hi-C matrix of a given experiment
    :param perc: percentage of cells with no count allowed to consider a column
       as valid.

    :returns: a dicitionary, which has as keys the index of the filtered out
       columns.
    """
    bads = {}
    size = len(matrx)
    cols = [size for i in xrange(size)]
    for k in matrx:
        cols[k / size] -= 1
    min_val = int(size * float(perc_zero) / 100)
    for i, col in enumerate(cols):
        if col > min_val:
            bads[i] = True
    if bads and not silent:
        stderr.write(('\nWARNING: removing columns having more than %s ' +
                      'zeroes:\n %s\n') % (
                         min_val, ' '.join(
                             ['%5s'%str(i + 1) + (''if (j + 1) % 20 else '\n')
                              for j, i in enumerate(sorted(bads.keys()))])))
    return bads


def hic_filtering_for_modelling(matrx, silent=False, perc_zero=90, auto=True,
                                draw_hist=False, savefig=None, diagonal=True):
    """
    Call filtering function, to remove artefactual columns in a given Hi-C
    matrix. This function will detect columns with very low interaction
    counts; and columns with NaN values (in this case NaN will be replaced
    by zero in the original Hi-C data matrix). Filtered out columns will be
    stored in the dictionary Experiment._zeros.
    
    :param matrx: Hi-C matrix of a given experiment
    :param False silent: does not warn for removed columns
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param True diagonal: remove row/columns with zero in the diagonal
    :param 90 perc_zero: maximum percentage of cells with no interactions
       allowed.
    :param True auto: if False, only filters based on the given percentage
       zeros

    :returns: the indexes of the columns not to be considered for the
       calculation of the z-score
    """
    bads = filter_by_zero_count(matrx, perc_zero, silent=silent)
    if auto:
        bads.update(filter_by_mean(matrx, draw_hist=draw_hist, silent=silent,
                                   savefig=savefig, bads=bads))
    # also removes rows or columns containing a NaN
    has_nans = False
    for i in xrange(len(matrx)):
        if matrx.get(i + i * len(matrx), 0) == 0 and diagonal:
            if not i in bads:
                bads[i] = None
        elif repr(sum([matrx.get(i + j * len(matrx), 0)
                       for j in xrange(len(matrx))])) == 'nan':
            has_nans = True
            if not i in bads:
                bads[i] = None
    return bads, has_nans
