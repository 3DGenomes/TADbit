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
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests

from skmisc.loess import loess
from serutils.stats.correlate import get_confidence, fit_with_uncertainty, latex_formula


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
    axe = inset_axes(ax, bbox_transform=ax.transAxes, height="100%",
                     width="100%", loc=3, bbox_to_anchor=(0, 0.1, 0.87, 0.87))
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
        leg = axe.legend(dots, [
            '%.1f%% of the bins more different than %.1f%% of null model)' %
            (len(signx) * 100. / total_len, cut * 100.)],
                         # bbox_to_anchor=(0.9, 1.03),
                         frameon=False)
    return axe


def nice_ba_plot(x, y, unadj_pN, sigmaN, sigmaR, pred, cond1, cond2,
                 alpha=0.75, ax=None):
    """
    Modified Bland-Altman plot (no logs) with prediction bands based on ordinary
     least square on comparison and null model.

    :param x: values average of eigenvectors
    :param y: values from difference of eigenvector
    :param x_cor: values average of eigenvectors (null model)
    :param y_cor: values difference of eigenvectors (null model)
    :param cond1: name of first experiment (corresponding to first eigenvector)
    :param cond2: name of second experiment (corresponding to second
       eigenvector)
    :param float('nan') alpha: used for LOESS fitting (only for labels)
    :param 0 df: number of degrees of freedom used in the fit
    :param None ax: matplotlib Axe object
    """
    plt.axis('off')
    vmin, vmax = -25, 0
    axe = inset_axes(ax, bbox_transform=ax.transAxes, height="100%",
                    width="100%", loc=3, bbox_to_anchor=(0, 0.1, 0.87, 0.87))
    axdw = inset_axes(ax, width="3%", height="45%", loc=3,
                    bbox_to_anchor=(0.89, 0.11, 0.98, 0.87),
                    bbox_transform=ax.transAxes, borderpad=0)
    axup = inset_axes(ax, width="3%", height="45%",
                        axes_kwargs={'facecolor': 'red'},
                        bbox_to_anchor=(0.89, 0.11, 0.98, 0.87),
                        loc=2, bbox_transform=ax.transAxes, borderpad=0)

    scup = axe.scatter(x[(y - pred > 0) & (unadj_pN < 0.05)],
                    y[(y - pred > 0) & (unadj_pN < 0.05)],
                    c=np.log(unadj_pN[(y - pred > 0) & (unadj_pN < 0.05)]),
                    s=50, alpha=1, facecolors='none', lw=0.5,
                    cmap='Oranges_r', marker='o', vmin=vmin, vmax=vmax)
    scdw = axe.scatter(x[(y - pred < 0) & (unadj_pN<0.05)],
                    y[(y - pred < 0) & (unadj_pN<0.05)],
                    c=np.log(unadj_pN[(y - pred < 0) & (unadj_pN < 0.05)]),
                    s=50, alpha=1, facecolors='none', lw=0.5,
                    cmap='Blues_r', marker='o', vmin=vmin, vmax=vmax)

    axe.plot(x, y, 'k.', alpha=0.6, ms=0.5)

    lines = []
    lines.extend(axe.plot(x, pred, 'k-',
                          label='LOESS fit ($\\alpha=%.1f$)' % (alpha)))
    lines.extend(axe.plot(x, pred - 2 * sigmaN, ':', color='grey', lw=1.5,
                        label='LOESS fit $\\pm2\\sigma$ from Null model'))
    axe.plot(x, pred + 2 * sigmaN, ':', color='grey', lw=1.5)
    lines.extend(axe.plot(x, pred - 2 * sigmaR, ':', color='black', lw=1.5,
                        label='LOESS fit $\\pm2\\sigma$ from comparison'))
    axe.plot(x, pred + 2 * sigmaR, ':', color='black'  , lw=1.5)
    axe.set_xlabel('Average of %s and %s Eigenvectors' % (cond1, cond2))
    axe.set_ylabel("Difference of Eigenvectors (%s - %s)" % (cond1, cond2))


    tickvals = np.linspace(0, vmin, 6)

    caxup = plt.colorbar(scup, cax=axup)
    caxup.set_ticks(tickvals)
    caxup.ax.hlines(1 - np.log(0.05) / vmin, 0, 1,
                    color='black', linestyle=':')
    caxup.ax.hlines(1 - np.log((sigmaR / sigmaN) * 2 *0.05) / vmin, 0, 1,
                    color='grey', linestyle=':')
    caxup.set_label('Log p-value student outlier test\n(for negative differences)')

    caxdw = plt.colorbar(scdw, cax=axdw)
    caxdw.set_ticks(tickvals)
    caxdw.ax.hlines(1 - np.log(0.05) / vmin, 0, 1,
                    color='black', linestyle=':')
    caxdw.ax.hlines(1 - np.log((sigmaR / sigmaN) * 2 * 0.05) / vmin, 0, 1,
                    color='grey', linestyle=':')
    caxdw.set_label('Log p-value student outlier test\n(for positive differences)')

    axe.legend(lines, [l.get_label() for l in lines], frameon=False,
               fontsize=9)

    return axe

def compare_AB(ev1, ev2, axe=None, xlabel='', ylabel='', color_ab=False,
               kde=True, use_odr=True, EV_range=(-1, 1), mid_point=0):
    if not axe:
        axe = plt.subplot(111)
    dots = []
    num_dots = []
    if color_ab:
        try:
            a, b = zip(*[(ev1[i], ev2[i]) for i in xrange(len(ev1)) if ev1[i] > mid_point and ev2[i] > mid_point])
            dots.extend(plt.plot(a, b, 'r.', alpha=0.1, label='Bins always in A'))
            num_dots.append(len(a))
        except ValueError:
            pass
        try:
            a, b = zip(*[(ev1[i], ev2[i]) for i in xrange(len(ev1)) if ev1[i] < mid_point and ev2[i] < mid_point])
            dots.extend(plt.plot(a, b, 'b.', alpha=0.1, label='Bins always in B'))
            num_dots.append(len(a))
        except ValueError:
            pass
        try:
            a, b = zip(*[(ev1[i], ev2[i]) for i in xrange(len(ev1))
                         if (ev1[i] < mid_point and ev2[i] > mid_point)
                         or (ev1[i] > mid_point and ev2[i] < mid_point)])
            dots.extend(plt.plot(a, b, '.', color='grey', alpha=0.1, label='Bins switching'))
            num_dots.append(len(a))
        except ValueError:
            pass
    else:
        dots.extend(plt.plot(ev1, ev2, '.', color='grey', alpha=0.1, label='Bins'))

    r, p  = st.pearsonr(ev1, ev2)
    plt.xlim(EV_range)
    plt.ylim(EV_range)
    plt.xticks([EV_range[0], mid_point, EV_range[1]])
    plt.yticks([EV_range[0], mid_point, EV_range[1]])
    plt.axhline(mid_point, color='lightgrey', alpha=0.3)
    plt.axvline(mid_point, color='lightgrey', alpha=0.3)
    p_x, p_y, z, confs, preds, r2 = fit_with_uncertainty(
        ev1, ev2, x_range=(-2, 2), use_odr=use_odr)
    formula = latex_formula("A*x+B", z)
    # confs, preds, p_x, p_y, z, r2, formula
    corr_color = "#7f7f7f"
    fit_line = plt.plot(p_x, p_y,color=corr_color, alpha=0.7, lw=2, label='Regression line')
    # plot confidence limits
    plt.fill_between(p_x, p_y - preds, p_y + preds, color=corr_color, alpha=0.15)
    plt.plot(p_x, p_y - preds, color=corr_color, alpha=0.15)
    plt.plot(p_x, p_y + preds, color=corr_color, alpha=0.25)

    p1 = Rectangle((0, 0), 1, 1, fc=corr_color, alpha=.2)
    num_dots = [100 * float(n) / sum(num_dots) for n in num_dots]
    dot_labels = (['Bin stays A: %.0f%%' % (num_dots[0]),
                   'Bin stays B: %.0f%%' % (num_dots[1]),
                   'Bin changes: %.0f%%' % (num_dots[2])]
                  if color_ab else ['Bins'])
    leg = plt.legend(fit_line + [p1] + dots,
                     ['''$y = %s$ (Pearson %.2f)''' % (formula, r),
                      '95% Prediction band'] + dot_labels,
                     frameon=False, loc=2)
    # a little bit of cheating here to see the dots
    for l in leg.get_lines():
        if not 'Bins' in l.get_label():
            continue
        l.set_alpha(0.6)
        x0, x1 = l.get_xdata()
        l.set_xdata([(x0 + x1) / 2])
        y0, y1 = l.get_ydata()
        l.set_ydata([(y0 + y1) / 2])
        l.set_marker('o')
        l.set_markersize(l.get_markersize() * 1.25)

    # plt.text(-0.7, 0.9, '$y = %s$' % formula)
    # plt.text(-0.7, 0.75, "Pearson=%.3f" % r)

    divider = make_axes_locatable(axe)
    # now determine nice limits by hand:
    axHistx = divider.append_axes("top", size=0.4, pad=0.025, sharex=axe)
    axHisty = divider.append_axes("right", size=0.4, pad=0.025, sharey=axe)

    if kde:
        kde1 = st.gaussian_kde(ev1)
        kde2 = st.gaussian_kde(ev2)

        xx1 = np.linspace(-1.1, 1.1, 100)
        axHistx.fill_between(xx1, kde1(xx1), color='LightGray')
        med1 = np.median(ev1)
        axHistx.vlines(med1, 0, kde1(med1), linestyle='-', color='grey')
        xx2 = np.linspace(-1.1, 1.1, 100)
        axHisty.fill_betweenx(xx2, kde2(xx2), color='LightGray')
        med2 = np.median(ev2)
        axHisty.hlines(med2, 0, kde2(med2), linestyle='-', color='grey')
    else:
        _ = axHistx.hist(ev1, bins=100)
        _ = axHisty.hist(ev2, bins=100, orientation='horizontal')

    axHisty.axison = False
    axHistx.axison = False
    axe.set_xlabel(xlabel)
    axe.set_ylabel(ylabel)
    axe.set_xlim(EV_range)
    axe.set_ylim(EV_range)
    return axHistx


def ba_plot(x, y, pred, cond1, cond2, alpha=float('nan'), df=0, ax=None):
    """
    Modified Bland-Altman plot (no logs)

    :param x: values average of eigenvectors
    :param y: values from difference of eigenvector
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

    # ~normalize
    ev1 = ev1 / np.std(ev1) / 3
    ev2 = ev2 / np.std(ev2) / 3

    # plot
    axes = []
    if plot == 'all':
        _ = plt.figure(figsize=(18, 18))
    elif plot != 'none':
        _ = plt.figure(figsize=(10, 10))
        axe = plt.subplot(111)
    if plot in ['all', 'correlation']:
        if plot == 'all':
            axe = plt.subplot(2, 2, 1)
        axes.append(axe)
        axe = compare_AB(ev1, ev2, axe=axe, xlabel='Eigenvector of ' + cond1,
                   ylabel='Eigenvector of ' + cond2, color_ab=True)
        axe.set_title('Correlation between EigenVectors (%s vs %s)' % (
            cond1, cond2))

    # Definition of null model
    print 'Defining Null model'
    ev3 = np.array(list(ev1[1:]) + list(ev2[1:]))
    ev4 = np.array([ev1[v] for v in xrange(len(ev1) - 1)] +
                   [ev2[v] for v in xrange(len(ev2) - 1)])

    if plot == 'all':
        axes.append(plt.subplot(2, 2, 2))
        axe = compare_AB(
            ev3, ev4, axe=axes[-1],
            xlabel='Eigenvector from %s and %s ($n$)' % (cond1, cond2),
            ylabel='Eigenvector from %s and %s ($n+1$)' % (cond1, cond2),
            color_ab=True)
        axe.set_title((
            'Correlation of EigenVectors (Null model)\n'
            '{0} bins $n$ vs $n+1$ and {1} '
            'bins $n$ vs $n+1$').format(cond1, cond2))

    ##########################################################################
    # Normalization

    # Z-scorish
    zev1 = (ev1 - np.mean(ev1)) / np.std(ev1) / 3
    zev2 = (ev2 - np.mean(ev2)) / np.std(ev2) / 3
    # prepare data for MA plot
    x = (zev1 + zev2) / 2
    y = (zev1 - zev2)
    # sort
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    ids = ids[idx]

    # for null model:
    zev3 = np.array(list(zev1[1:]) + list(zev2[1:]))
    zev4 = np.array([zev1[v] for v in xrange(len(zev1) - 1)] +
                   [zev2[v] for v in xrange(len(zev2) - 1)])
    x_cor = (zev3 + zev4) / 2
    y_cor = (zev3 - zev4)
    idx_cor = np.argsort(x_cor)
    x_cor = x_cor[idx_cor]
    y_cor = y_cor[idx_cor]

    # LOESS
    if norm == 'loess':
        print 'Computing LOESS on observed data'
        lo = loess(x, y, span=alpha, weight=None)
        lo.fit()
        pred = lo.outputs.fitted_values.copy()
        df = lo.outputs.enp
    else:
        pred = np.zeros(len(y))
        df = 0

    # LOESS on Null model
    if norm == 'loess':
        print 'Computing LOESS on Null model'
        lo = loess(x_cor, y_cor, span=alpha, weight=None)
        lo.fit()
        pred_cor = lo.outputs.fitted_values.copy()
        df = lo.outputs.enp
    else:
        pred_cor = np.zeros(len(y_cor))
        df = 0

    ##########################################################################
    # ordinary least square regression
    print 'Perform OLS regression and outlier test'

    modelR = OLS(y - pred, x, )
    resultR = modelR.fit()

    modelN = OLS(y_cor - pred_cor, x_cor)
    resultN = modelN.fit()

    sigmaN = np.sqrt(resultN.mse_resid)

    inflR = resultR.get_influence
    hiiR = inflR().hat_matrix_diag
    sigmaR = np.sqrt(resultR.mse_resid)
    residR = resultR.resid / sigmaN / np.sqrt(1 - hiiR)
    dfR = modelR.df_resid - 1

    unadj_pR = st.t.sf(np.abs(residR), dfR) * 2
    adj_pR = multipletests(unadj_pR, alpha=0.05, method='bonferroni')

    unadj_pN = st.t.sf(np.abs(residR), dfR) * 2
    adj_pN = multipletests(unadj_pN, alpha=0.05, method='bonferroni')

    if plot in ['all', 'difference']:
        if plot == 'all':
            axe = plt.subplot(2, 2, 3)
        axe.set_title(('Bland-Altman plot of EigenVectors (%s vs %s)\n'
                       'with prediction bands based on null model') % (
                           cond1, cond2))
        axes.append(nice_ba_plot(x, y, unadj_pN, sigmaN, sigmaR, pred,
                                 cond1, cond2, alpha=alpha, ax=axe))

    ##########################################################################
    print 'Perform the kernel density estimate for null model'
    y -= pred
    y_cor -= pred_cor
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

    ev1 = ev1[idx]
    ev2 = ev2[idx]
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
        result[ids[i][0], ids[i][1]] = (ev1[i], ev2[i], pv, y[i],
                                        unadj_pR[i], unadj_pN[i])

    if plot in ['all', 'density']:
        if plot == 'all':
            axe = plt.subplot(2, 2, 4)
        plt.title(('LOESS normalized BA density plot of EigenVectors\n'
                   '({0} vs {1} plotted over null model)').format(
                       cond1, cond2))
        axes.append(nice_contour_plot(
            xx, yy, f, f_cor, cond1, cond2, ax=axe, total_len=len(x), cut=cut,
            signx=signx, signy=signy, t_contours=t_contours, steps=steps))

    if plot in ['all', 'correlation']:
        xlim = (min(ev1.min(), ev3.min()), max(ev1.max(), ev3.max()))
        maxval = max(abs(xlim[0]), abs(xlim[1]))
        xlim = (-maxval, maxval)
        ylim = (min(ev2.min(), ev4.min()), max(ev2.max(), ev4.max()))
        maxval = max(abs(ylim[0]), abs(ylim[1]))
        ylim = (-maxval, maxval)
        for axe in axes[:2]:
            axe.set_xlim(xlim)
            axe.set_ylim(ylim)
    if plot in ['all', 'density', 'difference']:
        xlim = (min(x.min(), x_cor.min()), max(x.max(), x_cor.max()))
        ylim = (min(y.min(), y_cor.min()), max(y.max(), y_cor.max()))
        for axe in (axes[2:] if plot == 'all' else axes):
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

    result = get_significant_ev(ev1, ev2, names, opts.name1, opts.name2,
                                alpha=opts.alpha, confidence=opts.conf,
                                plot=opts.plot, norm=opts.norm,
                                kernel_density=opts.kernel_density)
    os.system('mkdir -p %s' % opts.outdir)
    if opts.plot != 'none':
        plt.savefig(os.path.join(opts.outdir,
                                 'differential_analysis_plot_%s_%s-%s.%s' % (
                                     opts.plot, opts.name1, opts.name2, opts.format)),
                    format=opts.format)

    out = open(os.path.join(opts.outdir,
                            'comparison_%s-%s.tsv' % (opts.name1,
                                                      opts.name2)), 'w')
    out.write('# Chromosome\tposition\t'
              'EV1/(3*std)\tEV2/(3*std)\t'
              'probability KDE vs null\t'
              'LOESS normalized differences\t'
              'student outlier test p-value OLS\t'
              'student outlier test p-value OLS (sigma from null model)\n')
    out.write('\n'.join(
        '%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f' % (a, b, i, j, pv, diff, pr, pn)
        for (a, b), (i, j, pv, diff, pr, pn) in result.iteritems()) + '\n')
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
                        default='all', choices=['all', 'none', 'density',
                                                'correlation', 'difference'],
                        help='''[%(default)s] plot all  comparisons, or only
                        one (density or difference), or none.''')
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
