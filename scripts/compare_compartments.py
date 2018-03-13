#!/usr/bin/env python
"""
Compare eigenvectors of A/B compartments between two experiments
"""

import os
from argparse import ArgumentParser

from scipy import interpolate
import scipy.stats as st

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from skmisc.loess import loess
from serutils.stats.correlate import get_confidence


def nice_contour_plot(xx, yy, f, f_cor, cond1, cond2, ax=None, t_contours=None,
                      steps=None, signx=None, signy=None, total_len=None,
                      cut=0.95):
    """
    Modified Bland-Altman density plot (no logs)

    :param xx: values from one eigenvector
    :param yy: values from other eigenvector
    :param f: kernel density function
    :param f_corr: kernel density function of null model
    :param cond1: name of first experiment (corresponding to first eigenvector)
    :param cond2: name of second experiment (corresponding to second
       eigenvector)
    :param None t_contours: steps at which to draw line in contour plot
    :param None steps: labels for steps at which to draw line in contour plot
    :param None signx: X coordinates of significant differences
    :param None signy: Y coordinates of significant differences
    :param None total_len: total number of bins
    :param 0.95 cut: cutoff to difine significant differences
    :param None ax: matplotlib Axe object
    """
    plt.axis('off')
    axe = inset_axes(ax,
                     width="100%",
                     height="100%",
                     loc=3,
                     bbox_transform=ax.transAxes,
                     bbox_to_anchor=(0.0, 0.1, 0.87, 0.87))
    if signx:
        dots = axe.plot(signx, signy, 'ko', mfc='none', ms=2.5, mew=0.5)
    # Contourf plot
    # plt.contour(z.T, t_contours, extent=[-3,3,-3,3])

    cfset = axe.contourf(xx, yy, f, np.linspace(0.0000001, f.max(), 200),
                         vmin=0.00001, locator=ticker.LogLocator(),
                         scale='log', cmap='Reds')
    axdown = inset_axes(ax, width="3%", height="45%", loc=3,
                        bbox_to_anchor=(0.89, 0.11, 0.98, 0.87),
                        bbox_transform=ax.transAxes, borderpad=0)
    cfbar = plt.colorbar(cfset, cax=axdown)
    cfbar.set_ticks(np.linspace(0.0000001, f.max(), 5))
    cfbar.set_ticks([])
    cfbar.set_ticklabels([])
    cfbar.ax.set_ylabel('Density of comparison obsevartions')
    if t_contours is not None:
        c_set = axe.contour(xx, yy, f_cor, t_contours, cmap='Greys_r',
                            locator=ticker.LogLocator())

        axup = inset_axes(ax, width="3%", height="45%",
                          axes_kwargs={'facecolor': 'red'},
                          bbox_to_anchor=(0.89, 0.11, 0.98, 0.87),
                          loc=2, bbox_transform=ax.transAxes, borderpad=0)
        c_bar = plt.colorbar(c_set, cax=axup, boundaries=(-0.1, 1.1))
        # cbar.set_ticks(linspace(0.0000001, f_cor.max(), 5))
        c_bar.set_ticklabels(['%d%%' % v for v in [99, 95, 90, 75, 50, 25]])
        lines = c_bar.ax.get_children()[0]
        lines.set_linewidths(7)
        c_bar.ax.set_ylabel('Percentages from Null model')
        c_bar.patch.set_facecolor('#FDD5C4')

        fmt = {}
        for l, s in zip(c_set.levels, steps):
            fmt[l] = str(int(s * 100)) + '%'
        axup.clabel(c_set, inline=1, fontsize=9, fmt=fmt, colors='k')

    # Label plot
    axe.set_ylabel("LOESS normalized difference of Eigenvectors (%s - %s)" % (
        cond1, cond2))
    axe.set_xlabel('Average of %s and %s Eigenvectors' % (cond1, cond2))
    if signx:
        ax.legend(dots, [
            '%.1f%% of the bins more different than %.1f%% of null model)' %
            (len(signx) * 100. / total_len, cut * 100.)],
                  bbox_to_anchor=(0.9, 1.03),
                  frameon=False)
    return axe


def ba_plot(x, y, pred, cond1, cond2, alpha=float('nan'), df=0, ax=None):
    """
    Modified Bland-Altman plot (no logs)

    :param x: values from one eigenvector
    :param y: values from other eigenvector
    :param cond1: name of first experiment (corresponding to first eigenvector)
    :param cond2: name of second experiment (corresponding to second
       eigenvector)
    :param float('nan') alpha: used for LOESS fitting (only for labels)
    :param 0 df: number of degrees of freedom used in the fit
    :param None ax: matplotlib Axe object
    """
    dots = ax.plot(x, y, 'k.', alpha=0.8, ms=0.5)
    fit_line = ax.plot(x, pred, label="loess (d = 1, f = %.2f)" % alpha,
                       color='red')
    confs, preds, _ = get_confidence(x, y, pred, x, df, conf=0.95)
    ax.fill_between(x, pred - preds, pred + preds, alpha=0.8, color='#FDD5C4')
    ax.fill_between(x, pred - confs, pred + confs, alpha=0.8, color='orange')
    ax.set_ylabel("Difference of Eigenvectors (%s - %s)" % (cond1, cond2))
    ax.set_xlabel('Average of %s and %s Eigenvectors' % (cond1, cond2))
    p1 = plt.Rectangle((0, 0), 1, 1, fc="#FDD5C4", alpha=.8)
    p2 = plt.Rectangle((0, 0), 1, 1, fc="orange", alpha=.8)
    ax.legend(dots + fit_line + [p1, p2],
              ['%s-%s pair' % (cond1, cond2),
               '''LOESS fit ($\\alpha=%.2f$)''' % (alpha),
               '95% Confidence band', '95% Prediction band'],
              loc='upper right', frameon=False, bbox_to_anchor=[1, 1])


def get_significant_ev(ev1, ev2, ids, cond1, cond2, norm='loess', plot='all',
                       confidence=0.95, alpha=0.75, kernel_density=200):
    """
    Compare two eigenvectors (from two conditions), using as null, or
       background, model the differences between each pair of neighbor bins.
       EigenVectors are Z-score normalized and their difference is LOESS
       normalized

    :param ev1: list of values in eigenvectors from first condition
    :param ev2: list of values in eigenvectors from second condition
    :param ids: list of names (e.g. [('chr1', 1294), ...]) corresponding to
       values in each eigenvector
    :param cond1: name of first experiment (corresponding to first eigenvector)
    :param cond2: name of second experiment (corresponding to second
       eigenvector)
    :param loess norm: normalization to perform on the difference of
       eigenvectors (options are loess or none)
    :param 0.75 alpha: smoothing parameter for LOESS normalization (0 pass
       through all points, 1, straight line)
    :param 200 kernel_density: density of the matrix for Gaussian kernel
       density estimation
    :param all plot: plots to be generate. If 'all' 6 plots  will be generated
       three for differences between conditions (before, after LOESS
       normalization, and density map), the same two for null mode. If 'final'
       only a single plot with density maps of observed data and null model. If
       'none', no plot will be generated.
    :param 0.95 confidence: confidence level for definition of bins with
       significantly different compartments

    :returns: a dictionary with, as keys, the input ids, and as values, a tuple
       with two value: 1- the probabilities (or Cumulative Densities) of each
       compared pair of eigenvector values inside the Gaussian kernel density
       of the null model, and 2- the LOESS normalized difference between
       eigenvectors values
    """
    print 'Getting EigenVectors'
    ev1 = np.array(ev1)
    ev2 = np.array(ev2)
    ids = np.array(ids)
    # prepare data for MA plot
    x = (ev1 + ev2) / 2
    y = (ev1 - ev2)
    # sort
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    ids = ids[idx]

    # LOESS
    if norm == 'loess':
        print 'Computing LOESS on observed data'
        lo = loess(x, y, span=alpha, weight=None)
        lo.fit()
        pred = lo.outputs.fitted_values
        df = lo.outputs.enp
    else:
        pred = np.zeros(len(y))
        df = 0
    # plot
    axes = []
    if plot == 'all':
        _ = plt.figure(figsize=(18, 27))
    elif plot == 'final':
        _ = plt.figure(figsize=(10, 10))
    if plot == 'all':
        axes.append(plt.subplot(3, 2, 2))
        axes[-1].set_title('Bland-Altman plot of EigenVectors (%s vs %s)' % (
            cond1, cond2))
        ba_plot(x, y, pred, cond1, cond2, alpha, df, axes[-1])

    # LOESS normalization
    y = y - pred
    print 'Computing LOESS on normalized observed data'
    if norm == 'loess':
        lo = loess(x, y, span=alpha, weight=None)
        lo.fit()
        pred = lo.outputs.fitted_values
        df = lo.outputs.enp
    else:
        pred = np.zeros(len(y))
        df = 0
    # plot
    if plot == 'all':
        axes.append(plt.subplot(3, 2, 4))
        axes[-1].set_title(('LOESS normalized BA plot of EigenVectors '
                            '(%s vs %s)') % (cond1, cond2))
        ba_plot(x, y, pred, cond1, cond2, alpha, df, axes[-1])

    # Definition of null model
    print 'Defining Null model'
    ev3 = np.array(list(ev1[1:]) + list(ev2[1:]))
    ev4 = np.array([ev1[v] for v in xrange(len(ev1) - 1)] +
                   [ev2[v] for v in xrange(len(ev2) - 1)])
    x_cor = (ev3 + ev4) / 2
    y_cor = (ev3 - ev4)
    idx = np.argsort(x_cor)
    x_cor = x_cor[idx]
    y_cor = y_cor[idx]

    # LOESS on Null model
    if norm == 'loess':
        print 'Computing LOESS on Null model'
        lo = loess(x_cor, y_cor, span=alpha, weight=None)
        lo.fit()
        pred_cor = lo.outputs.fitted_values
        df = lo.outputs.enp
    else:
        pred_cor = np.zeros(len(y_cor))
        df = 0
    if plot == 'all':
        axes.append(plt.subplot(3, 2, 1))
        axes[-1].set_title((
            'Bland-Altman plot of EigenVectors (Null model)\n'
            '{0} even bins ($n$) vs odd ($n+1$) and {1} even '
            'bins ($n$) vs odd ($n+1$)\n').format(cond1, cond2))
        ba_plot(x_cor, y_cor, pred_cor, cond1, cond2, alpha, df, axes[-1])

    # LOESS normalization
    y_cor = y_cor - pred_cor
    print 'Computing LOESS on normalized observed data'
    if norm == 'loess':
        lo = loess(x_cor, y_cor, span=alpha, weight=None)
        lo.fit()
        pred_cor = lo.outputs.fitted_values
        df = lo.outputs.enp
    else:
        pred_cor = np.zeros(len(y_cor))
        df = 0
    # plot
    if plot == 'all':
        axes.append(plt.subplot(3, 2, 3))
        axes[-1].set_title(
            ('LOESS normalized BA plot of EigenVectors (Null model)\n'
             '{0} and {1} together, even bins ($n$) vs odd ($n+1$)\n').format(
                 cond1, cond2))
        ba_plot(x_cor, y_cor, pred_cor, cond1, cond2, alpha, df, axes[-1])

    print 'Perform the kernel density estimate for null model'
    # Perform the kernel density estimate for null model
    xmin = min(x_cor) * 1.5
    ymin = min(y_cor) * 1.5
    xmax = max(x_cor) * 1.5
    ymax = max(y_cor) * 1.5
    xx, yy = np.mgrid[xmin:xmax:complex(0, kernel_density),
                      ymin:ymax:complex(0, kernel_density)]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    z_cor = np.vstack([x_cor, y_cor])
    kernel_cor = st.gaussian_kde(z_cor)
    f_cor = np.reshape(kernel_cor(positions).T, xx.shape)
    f_cor_sum = f_cor.sum()
    f_cor /= f_cor_sum

    print 'Perform the kernel density estimate for comparison'
    # Perform the kernel density estimate for comparison
    z = np.vstack([x, y])
    kernel = st.gaussian_kde(z)
    f = np.reshape(kernel(positions).T, xx.shape)
    f /= f.sum()

    # define probability lines
    n = 10000
    t = np.linspace(0, f_cor.max(), n)

    # kernel probability for null model
    integral = ((f_cor >= t[:, None, None]) * f_cor).sum(axis=(1, 2))
    # function to get kernel density at a given CDF
    ff = interpolate.interp1d(integral, t)
    steps = [0.99, 0.95, 0.9, 0.75, 0.50, 0.25]
    t_contours = ff(np.array(steps))
    # function to get CDF at a given kernel density
    invff = interpolate.interp1d(t, integral)

    # significant bins in Null model
    cut = confidence

    # significant bins observed data
    print 'Computing significant changes in observed data'
    signx = []
    signy = []
    result = {}

    # get kernel density for each pair of eigenvectors
    pvals = kernel_cor.pdf((x, y)) / f_cor_sum

    for i, pv in enumerate(pvals):
        try:
            pv = invff(pv)  # convert PDF to CDF
        except ValueError:
            try:
                pv = 1.
            except ValueError:
                pv = 0.
            continue
        if pv > cut:
            signx.append(x[i])
            signy.append(y[i])
        result[ids[i][0], ids[i][1]] = pv, y[i]

    # plot Null model
    if plot == 'all':
        axe = plt.subplot(3, 2, 5)
        print 'Plotting densities'
        plt.title(('LOESS normalized BA density plot of EigenVectors (Null '
                   'model)\n{0} and {1} together, even bins ($n$) vs odd '
                   '($n+1$)\n').format(cond1, cond2))
        axes.append(nice_contour_plot(
            xx, yy, f_cor, f_cor, cond1, cond2, ax=axe,
            t_contours=t_contours, steps=steps))
        # plot Observed data
    if plot in ['all', 'final']:
        if plot == 'all':
            axe = plt.subplot(3, 2, 6)
        else:
            axe = plt.subplot(111)
        plt.title(('LOESS normalized BA density plot of EigenVectors\n'
                   '({0} vs {1} plotted over null model)\n').format(
                       cond1, cond2))
        axes.append(nice_contour_plot(
            xx, yy, f, f_cor, cond1, cond2, ax=axe, total_len=len(x), cut=cut,
            signx=signx, signy=signy, t_contours=t_contours, steps=steps))

        xlim = (min((x.min(), x_cor.min())), max((x.max(), x_cor.max())))
        ylim = (min((y.min(), y_cor.min())), max((y.max(), y_cor.max())))
        for axe in axes:
            axe.set_xlim(xlim)
            axe.set_ylim(ylim)
    return result


def main():
    """
    Main function for command line
    """
    opts = get_options()
    names = []
    ev1 = []
    ev2 = []
    for fnam in opts.fnam1:
        chrom = os.path.split(fnam)[-1].split('_')[0]
        fh = open(fnam)
        fh.next()
        ev1.extend(float(line.split()[0]) for line in fh)
        names.extend((chrom, i) for i in xrange(0, len(ev1) - len(names)))

    if len(names) != len(set(names)):
        raise Exception('ERROR: found duplicated entries.')

    for fnam in opts.fnam2:
        fh = open(fnam)
        chrom = os.path.split(fnam)[-1].split('_')[0]
        fh.next()
        ev2.extend(float(line.split()[0]) for line in fh)

    if len(ev1) != len(ev2):
        raise Exception('ERROR: number of bins differ between experiments.')

    # remove NaNs:
    for i in xrange(len(ev1) - 1, -1, -1):
        if not np.isfinite(ev1[i] + ev2[i]):
            ev1.pop(i)
            ev2.pop(i)
            names.pop(i)

    # Z-scorish
    ev1 = (ev1 - np.mean(ev1)) / np.std(ev1) / 3
    ev2 = (ev2 - np.mean(ev2)) / np.std(ev2) / 3

    result = get_significant_ev(ev1, ev2, names, opts.name1, opts.name2,
                                alpha=opts.alpha, confidence=opts.conf,
                                plot=opts.plot, norm=opts.norm,
                                kernel_density=opts.kernel_density)
    os.system('mkdir -p %s' % opts.outdir)
    plt.savefig(os.path.join(opts.outdir,
                             'differential_analysis_plot_%s-%s.%s' % (
                                 opts.name1, opts.name2, opts.format)),
                format=opts.format)

    out = open(os.path.join(opts.outdir,
                            'comparison_%s-%s.tsv' % (opts.name1,
                                                      opts.name2)), 'w')
    out.write('\n'.join('%s\t%s\t%f\t%f' % (a, b, pv, diff)
                        for (a, b), (pv, diff) in result.iteritems()) + '\n')
    out.close()


def get_options():
    """
    Options for command line
    """
    parser = ArgumentParser()

    parser.add_argument('--ev1', dest='fnam1', metavar='PATH', required=True,
                        default=False, nargs='+',
                        help='input list 1 of files with EigenVectors')
    parser.add_argument('--ev2', dest='fnam2', metavar='PATH', required=True,
                        default=False, nargs='+',
                        help='input list 2 of files with EigenVectors')
    parser.add_argument('--name1', dest='name1', metavar='STR',
                        default='Condition1',
                        help='name of input list 1 of EigenVectors')
    parser.add_argument('--name2', dest='name2', metavar='STR',
                        default='Condition2',
                        help='name of input list 2 of EigenVectors')
    parser.add_argument('-r', '--resolution', dest='reso', required=True,
                        type=int,
                        help='resolution at which compartments were called')
    parser.add_argument('-o', dest='outdir', metavar='PATH',
                        default='tmp/',
                        help='[%(default)s] directory to store results')
    parser.add_argument('--plot', dest='plot',
                        default='all', choices=['all', 'none', 'final'],
                        help='''[%(default)s] plot all stages of the
                        comparison, or only the final stage, are none.''')
    parser.add_argument('--format', dest='format', default='png',
                        help='''[%(default)s] plot format''')
    parser.add_argument('--norm', dest='norm',
                        default='loess', choices=['loess', 'none'],
                        help='''[%(default)s] flatten differences between
                        conditions.''')
    parser.add_argument('-k', '--kernel_density', dest='kernel_density',
                        default=100, type=int,
                        help='''[%(default)s] density of the matrix for
                        Gaussian kernel density estimation.''')
    parser.add_argument('-a', '--alpha', dest='alpha',  default=0.75,
                        type=float,
                        help='''[%(default)s] alpha (smoothing parameter) for
                        LOESS fitting (0 pass through all points, 1, straight
                        line).''')
    parser.add_argument('--conf', dest='conf',  default=0.95,
                        type=float,
                        help='''[%(default)s] confidence level for definition
                        of bins with significantly different compartments.''')
    opts = parser.parse_args()
    return opts


if __name__ == '__main__':
    exit(main())
