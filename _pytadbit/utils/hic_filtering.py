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


def filter_by_zero_count(matrx, perc_zero, min_count=None, silent=True):
    """
    :param matrx: Hi-C matrix of a given experiment
    :param perc: percentage of cells with no count allowed to consider a column
       as valid.
    :param None min_count: minimum number of reads mapped to a bin (recommended
       value could be 2500). If set this option overrides the perc_zero
       filtering... This option is slightly slower.

    :returns: a dicitionary, which has as keys the index of the filtered out
       columns.
    """
    bads = {}
    size = len(matrx)
    if min_count is None:
        cols = [size for i in xrange(size)]
        for k in matrx:
            cols[k / size] -= 1
        min_val = int(size * float(perc_zero) / 100)
    else:
        if matrx.symmetricized:
            min_count *= 2
        cols = [0 for i in xrange(size)]
        for k, v in matrx.iteritems(): # linear representation of the matrix
            cols[k / size] += v
        min_val = size - min_count
    if min_count is None:
        check = lambda x: x > min_val
    else:
        check = lambda x: x < min_count
    for i, col in enumerate(cols):
        if check(col):
            bads[i] = True
    if bads and not silent:
        if min_count is None:
            stderr.write(('\nWARNING: removing columns having more than %s ' +
                          'zeroes:\n %s\n') % (
                             min_val, ' '.join(
                                 ['%5s' % str(i + 1) + (''if (j + 1) % 20 else '\n')
                                  for j, i in enumerate(sorted(bads.keys()))])))
        else:
            stderr.write(('\nWARNING: removing columns having less than %s ' +
                          'counts:\n %s\n') % (
                             int(size - min_val), ' '.join(
                                 ['%5s' % str(i + 1) + (''if (j + 1) % 20 else '\n')
                                  for j, i in enumerate(sorted(bads.keys()))])))
    return bads


def hic_filtering_for_modelling(matrx, silent=False, perc_zero=90, auto=True,
                                min_count=None, draw_hist=False, savefig=None,
                                diagonal=True):
    """
    Call filtering function, to remove artifactual columns in a given Hi-C
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
    :param None min_count: minimum number of reads mapped to a bin (recommended
       value could be 2500). If set this option overrides the perc_zero
       filtering... This option is slightly slower.
    :param True auto: if False, only filters based on the given percentage
       zeros

    :returns: the indexes of the columns not to be considered for the
       calculation of the z-score
    """
    bads = filter_by_zero_count(matrx, perc_zero, min_count=min_count, silent=silent)
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


def _best_window_size(sorted_prc, size, beg, end, verbose=False):
    """
    Search for best window size.
    Between given begin and end percentiles of the distribution of cis interactions
    searches for a window size (number of bins) where all median values are between
    median * stddev and median * stddev of the global measure.

    :param sorted_prc: list of percentages of cis interactions by bins, sorted
       by the total interactions in the corresponding bins.
    :param size: total number of bins
    :param beg: starting position of the region with expected 'normal' behavior
       of the cis-percentage
    :param end: last position of the region with expected 'normal' behavior
       of the cis-percentage
    :param False verbose: print running information

    :returns: window size
    """
    if verbose:
        print ('      -> defining window in number of bins to average values of\n'
               '         percentage of cis interactions')
    nwins = min((1000, size / 10))
    if nwins < 100:
        warn('WARNING: matrix probably too small to automatically filter out bins\n')
    win_size = 0
    prevn = 0
    count = 0
    # iterate over possible window sizes (use logspace to gain some time)
    for n in np.logspace(1, 4, num=100):
        n = int(n)
        if n == prevn:
            continue
        prevn = n

        tmp_std = []
        tmp_med = []
        for k in xrange(int(size * beg), int(size * end),
                        (int(size * end) - int(size * beg)) / nwins):
            vals = sorted_prc[k:k + n]
            tmp_std.append(np.std(vals))
            tmp_med.append(np.median(vals))
        med_mid = np.median([tmp_med[i] for i in xrange(nwins)])
        results = [m - s < med_mid < m + s
                   for m, s in zip(tmp_med, tmp_std)]

        # if verbose:
        #     print '       -', n, med_mid, sum(results)
        if all(results):
            if not count:
                win_size = n
            count += 1
            if count == 10:
                break
        else:
            count = 0

    if verbose:
        print '        * first window size with stable median of cis-percentage: %d' % (win_size)
    return win_size


def filter_by_cis_percentage(cisprc, beg=0.3, end=0.8, sigma=2, verbose=False,
                             savefig=None):
    """
    Define artifactual columns with either too low or too high counts of
    interactions by compraing their percentage of cis interactions
    (inter-chromosomal).

    :param cisprc: dictionary with counts of cis-percentage by bin number.
       Values of the dictionary are tuple with,m as first element the number
       of cis interactions and as second element the total number of
       interactions.
    :param 0.3 beg: proportion of bins to be considered as possibly having low
       counts
    :param 0.8 end: proportion of bins to be considered as possibly having high
       counts
    :param 2 sigma: number of standard deviations used to define lower and upper
       ranges in the varaition of the percentage of cis interactions
    :param None sevefig: path to save image of the distribution of cis
       percentages and total counts by bin.

    :returns: dictionary of bins to be filtered out (with either too low or too
       high counts of interactions).
    """
    sorted_sum, indices = zip(*sorted((cisprc[i][1], i) for i in cisprc))

    sorted_prc = [float(cisprc[i][0]) / cisprc[i][1] for i in indices]

    size = len(indices)

    win_size = _best_window_size(sorted_prc, size, beg, end, verbose=verbose)

    # define confidance bands, compute median plus/minus one standard deviation
    errors_pos = []
    errors_neg = []
    for k in xrange(0, size, 1):
        vals = sorted_prc[k:k+win_size]
        std = np.std(vals)
        med = np.median(vals)
        errors_pos.append(med + std * sigma)
        errors_neg.append(med - std * sigma)

    # calculate median and variation of median plus/minus one standard deviation
    # for values between percentile 10 and 90 of the distribution of the
    # percentage of cis interactions
    #  - for median plus one standard deviation
    std_err_pos =  np.std   (errors_pos[int(size * beg):int(size * end)])
    med_err_pos =  np.median(errors_pos[int(size * beg):int(size * end)])
    #  - for median minus one standard deviation
    std_err_neg =  np.std   (errors_neg[int(size * beg):int(size * end)])
    med_err_neg =  np.median(errors_neg[int(size * beg):int(size * end)])

    # define cutoffs, values of cis percentage plus 1 stddev should be between
    # the general median +/- 2 stddev of the distribution of the cis percentage
    # plus 1 stddev. Same on the side of median cis percentage minus 1 stddev
    beg_pos = med_err_pos - std_err_pos * sigma
    end_pos = med_err_pos + std_err_pos * sigma
    beg_neg = med_err_neg - std_err_neg * sigma
    end_neg = med_err_neg + std_err_neg * sigma

    cutoffL = None
    passed = 0
    consecutive = 10
    for cutoffL, (p, n) in enumerate(zip(errors_pos, errors_neg)):
        # print '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f' % (beg_pos, p, end_pos, beg_neg, n, end_neg)
        if (beg_pos < p < end_pos) and (beg_neg < n < end_neg):
            if passed >= consecutive:
                break
            passed += 1
        else:
            passed = 0
    else:
        raise Exception('ERROR: left cutoff not found!!!')
    cutoffL -= consecutive  # rescale, we asked for XX consecutive

    # right
    cutoffR = None
    passed = 0
    for cutoffR, (p, n) in enumerate(zip(errors_pos, errors_neg)[::-1]):
        cutoffR = size - cutoffR
        # print '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f' % (beg_pos, p, end_pos, beg_neg, n, end_neg)
        if (beg_pos < p < end_pos) and (beg_neg < n < end_neg):
            if passed >= consecutive:
                break
            passed += 1
        else:
            passed = 0
    else:
        raise Exception('ERROR: right cutoff not found!!!')
    cutoffR += consecutive  # rescale, we asked for XX consecutive


    min_count = sorted_sum[int(cutoffL)]
    max_count = sorted_sum[int(cutoffR)]

    if verbose:
        print '        * Lower cutoff applied until bin number: %d' % (cutoffL)
        print '        * too few  interactions defined as less than %9d interactions' % (
            min_count)
        print '        * Upper cutoff applied until bin number: %d' % (cutoffR)
        print '        * too much interactions defined as more than %9d interactions' % (
            max_count)

    # plot
    if savefig:
        if verbose:
            print '      -> Making plot...'
        fig = plt.figure(figsize=(20,11))
        ax1 = fig.add_subplot(111)
        plt.subplots_adjust(left=0.25, bottom=0.2)
        line1 = plt.plot(range(size),
                         [float(cisprc[i][0]) / cisprc[i][1] for i in indices],
                         '.', color='grey', alpha=0.2,
                         label='cis interactions ratio by bin', zorder=1)
        line2 = plt.plot(range(0, size, 20), [sum(float(cisprc[j][0]) / cisprc[j][1]
                                                  for j in indices[k:k+win_size]
                                                  if j in cisprc) / win_size
                                              for k in xrange(0, size, 20)],
                         '.', color='k', alpha=0.3,
                         label='cis interactions ratio by %d bin' % win_size,
                         zorder=1)

        for k, (p, n) in enumerate(zip(errors_pos[::size / 100], errors_neg[::size / 100])):
            plt.vlines(k * (size / 100), (p + n) / 2, p, color='red', alpha=0.6)
            plt.vlines(k * (size / 100), n, (p + n) / 2, color='blue', alpha=0.6)
        plt.plot(range(0, size, size / 100), errors_neg[::size/100], 'b^', mec='blue', alpha=0.5)
        plt.plot(range(0, size, size / 100), errors_pos[::size/100], 'rv', mec='red', alpha=0.5)

        plt.fill_between([0, size], beg_pos, end_pos, color='red', alpha=0.3, zorder=2)
        plt.text(-size/15., (end_pos + beg_pos) / 2, 'Confidance band for\nupper stddev of median',
                 color='red', ha='right', va='center')
        plt.fill_between([0, size], beg_neg, end_neg, color='blue', alpha=0.3, zorder=2)
        plt.text(-size/15., (end_neg + beg_neg) / 2, 'Confidance band for\nlower stddev of median',
                 color='blue', ha='right', va='center')

        plt.ylim((0,1.1))
        plt.ylabel('Ratio of cis interactions ratio')
        plt.fill_betweenx([0, 1.1], cutoffL, cutoffR, color='green', alpha=0.2)
        plt.text((cutoffR + cutoffL) / 2, -0.1,
                 ('Kept bins, top and bottom deviations from median cis-ratio\n' +
                  'should be inside their respective confidance bands'),
                 ha='center', color='green')

        ax2   = fig.add_subplot(111, sharex=ax1, frameon=False)
        line3 = ax2.plot(sorted_sum, 'rx', alpha=0.4, label='Log sum of interactions by bin')
        ax2.set_yscale('log')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_ylabel('Log interaction counts')
        lns   = line1 + line2 + line3
        labs  = [l.get_label() for l in lns]
        ax2.legend(lns, labs, loc=0, bbox_to_anchor=(0, 0), frameon=False)

        plt.title('Keeping from %.2f%% to %.2f%%' % (100 * float(cutoffL) / size,
                                                     100 * float(cutoffR) / size))
        plt.xlim(0, size)

        tadbit_savefig(savefig)
        plt.close('all')


    badcol = {}
    countL = 0
    countZ = 0
    countU = 0
    for c in xrange(size):
        if cisprc.get(c, [0, 0])[1] < min_count:
            badcol[c] = cisprc.get(c, [0, 0])[1]
            countL += 1
            if not c in cisprc:
                countZ += 1
        elif cisprc[c][1] > max_count:  # don't need get here, already cought in previous condition
            badcol[c] = cisprc.get(c, [0, 0])[1]
            countU += 1
    print '     => %d BAD bins (%d/%d/%d null/low/high counts) of %d (%.1f%%)' % (
        len(badcol), countZ, countL, countU, size, float(len(badcol)) / size * 100)

    return badcol
