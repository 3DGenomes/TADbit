"""
18 Nov 2014
"""

from pytadbit                     import HiC_data
from pytadbit.utils.extraviews    import tadbit_savefig, setup_plot
from pytadbit.utils.tadmaths      import nozero_log_matrix as nozero_log
from pytadbit.utils.tadmaths      import right_double_mad as mad
from warnings                     import warn
from collections                  import OrderedDict
from pytadbit.parsers.hic_parser  import load_hic_data_from_reads
from pytadbit.utils.extraviews    import nicer
from pytadbit.utils.file_handling import mkdir
from scipy.stats                  import norm as sc_norm, skew, kurtosis
from scipy.stats                  import pearsonr, spearmanr, linregress
from numpy.linalg                 import eigh
import numpy as np

try:
    from matplotlib import rcParams
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:
    warn('matplotlib not found\n')


def hic_map(data, resolution=None, normalized=False, masked=None,
            by_chrom=False, savefig=None, show=False, savedata=None,
            focus=None, clim=None, cmap='jet', pdf=False, decay=True,
            perc=10, name=None, decay_resolution=None, **kwargs):
    """
    function to retrieve data from HiC-data object. Data can be stored as
    a square matrix, or drawn using matplotlib

    :param data: can be either a path to a file with pre-processed reads
       (filtered or not), or a Hi-C-data object
    :param None resolution: at which to bin the data (try having a dense matrix
       with < 10% of cells with zero interaction counts). Note: not necessary
       if a hic_data object is passed as 'data'.
    :param False normalized: used normalized data, based on precalculated biases
    :param masked: a list of columns to be removed. Usually because to few
       interactions
    :param False by_chrom: data can be stored in a partitioned way. This
       parameter can take the values of:
        * 'intra': one output per each chromosome will be created
        * 'inter': one output per each possible pair of chromosome will be
           created
        * 'all'  : both of the above outputs
    :param None savefig: path where to store the output images. Note that, if
       the by_chrom option is used, then savefig will be the name of the
       directory containing the output files.
    :param None savedata: path where to store the output matrices. Note that, if
       the by_chrom option is used, then savefig will be the name of the
       directory containing the output files.
    :param None focus: can be either two number (i.e.: (1, 100)) specifying the
       start and end position of the sub-matrix to display (start and end, along
       the diagonal of the original matrix); or directly a chromosome name; or
       two chromosome names (i.e.: focus=('chr2, chrX')), in order to store the
       data corresponding to inter chromosomal interactions between these two
       chromosomes
    :param True decay: plot the correlation between genomic distance and
       interactions (usually a decay).
    :param False force_image: force to generate an image even if resolution is
       crazy...
    :param None clim: cutoff for the upper and lower bound in the coloring scale
       of the heatmap
    :param False pdf: when using the bny_chrom option, to specify the format of
       the stored images
    :param Reds cmap: color map to be used for the heatmap
    :param None decay_resolution: chromatin fragment size to consider when
       calculating decay of the number of interactions with genomic distance.
       Default is equal to resolution of the matrix.
    """
    if isinstance(data, str):
        data = load_hic_data_from_reads(data, resolution=resolution, **kwargs)
        if not kwargs.get('get_sections', True) and decay:
            warn('WARNING: not decay not available when get_sections is off.')
            decay = False
    hic_data = data
    resolution = data.resolution
    if not decay_resolution:
        decay_resolution = resolution
    if hic_data.bads and not masked:
        masked = hic_data.bads
    # save and draw the data
    if by_chrom:
        if focus:
            raise Exception('Incompatible options focus and by_chrom\n')
        if savedata:
            mkdir(savedata)
        if savefig:
            mkdir(savefig)
        for i, crm1 in enumerate(hic_data.chromosomes):
            for crm2 in hic_data.chromosomes.keys()[i:]:
                if by_chrom == 'intra' and crm1 != crm2:
                    continue
                if by_chrom == 'inter' and crm1 == crm2:
                    continue
                try:
                    subdata = hic_data.get_matrix(focus=(crm1, crm2), normalized=normalized)
                    start1, _ = hic_data.section_pos[crm1]
                    start2, _ = hic_data.section_pos[crm2]
                    masked1 = {}
                    masked2 = {}
                    if focus and hic_data.bads:
                        # rescale masked
                        masked1 = dict([(m - start1, hic_data.bads[m])
                                        for m in hic_data.bads])
                        masked2 = dict([(m - start2, hic_data.bads[m])
                                        for m in hic_data.bads])
                    if masked1 or masked2:
                        for i in xrange(len(subdata)):
                            if i in masked1:
                                subdata[i] = [float('nan')
                                              for j in xrange(len(subdata))]
                            for j in xrange(len(subdata)):
                                if j in masked2:
                                    subdata[i][j] = float('nan')
                    if savedata:
                        hic_data.write_matrix('%s/%s.mat' % (
                            savedata, '_'.join(set((crm1, crm2)))),
                                              focus=(crm1, crm2),
                                              normalized=normalized)
                    if show or savefig:
                        if (len(subdata) > 10000
                            and not kwargs.get('force_image', False)):
                            warn('WARNING: Matrix image not created, more than '
                                 '10000 rows, use a lower resolution to create images')
                            continue
                        draw_map(subdata, 
                                 OrderedDict([(k, hic_data.chromosomes[k])
                                              for k in hic_data.chromosomes.keys()
                                              if k in [crm1, crm2]]),
                                 hic_data.section_pos,
                                 '%s/%s.%s' % (savefig,
                                               '_'.join(set((crm1, crm2))),
                                               'pdf' if pdf else 'png'),
                                 show, one=True, clim=clim, cmap=cmap,
                                 decay_resolution=decay_resolution, perc=perc,
                                 name=name, cistrans=float('NaN'))
                except ValueError, e:
                    print 'Value ERROR: problem with chromosome %s' % crm1
                    print str(e)
                except IndexError, e:
                    print 'Index ERROR: problem with chromosome %s' % crm1
                    print str(e)
    else:
        if savedata:
            hic_data.write_matrix(savedata, focus=focus,
                                  normalized=normalized)
        if show or savefig:
            subdata = hic_data.get_matrix(focus=focus, normalized=normalized)
            if (len(subdata) > 10000 and not kwargs.get('force_image', False)):
                warn('WARNING: Matrix image not created, more than '
                     '10000 rows, use a lower resolution to create images')
                return
            start1 = hic_data._focus_coords(focus)[0]
            if focus and masked:
                # rescale masked
                masked = dict([(m - start1, masked[m]) for m in masked])
            if masked:
                for i in xrange(len(subdata)):
                    if i in masked:
                        subdata[i] = [float('nan')
                                      for j in xrange(len(subdata))]
                    for j in xrange(len(subdata)):
                        if j in masked:
                            subdata[i][j] = float('nan')
            draw_map(subdata,
                     {} if focus else hic_data.chromosomes,
                     hic_data.section_pos, savefig, show,
                     one = True if focus else False, decay=decay,
                     clim=clim, cmap=cmap, decay_resolution=decay_resolution,
                     perc=perc, normalized=normalized,
                     max_diff=kwargs.get('max_diff', None),
                     name=name, cistrans=float('NaN') if focus else
                     hic_data.cis_trans_ratio(normalized,
                                              kwargs.get('exclude', None),
                                              kwargs.get('diagonal', True),
                                              kwargs.get('equals', None)))


def draw_map(data, genome_seq, cumcs, savefig, show, one=False, clim=None,
             cmap='jet', decay=False, perc=10, name=None, cistrans=None,
             decay_resolution=10000, normalized=False, max_diff=None):
    _ = plt.figure(figsize=(15.,12.5))
    if not max_diff:
        max_diff = len(data)
    ax1 = plt.axes([0.34, 0.08, 0.6, 0.7205])
    ax2 = plt.axes([0.07, 0.65, 0.21, 0.15])
    if decay:
        ax3 = plt.axes([0.07, 0.42, 0.21, 0.15])
        plot_distance_vs_interactions(data, genome_seq=genome_seq, axe=ax3,
                                      resolution=decay_resolution,
                                      max_diff=max_diff, normalized=normalized)
    ax4 = plt.axes([0.34, 0.805, 0.6, 0.04], sharex=ax1)
    ax5 = plt.axes([0.34, 0.845, 0.6, 0.04], sharex=ax1)
    ax6 = plt.axes([0.34, 0.885, 0.6, 0.04], sharex=ax1)
    try:
        minoridata   = np.nanmin(data)
        maxoridata   = np.nanmax(data)
    except AttributeError:
        vals = [i for d in data for i in d if not np.isnan(i)]
        minoridata   = np.min(vals)
        maxoridata   = np.max(vals)
    totaloridata = np.nansum([data[i][j] for i in xrange(len(data))
                              for j in xrange(i, len(data[i]))]) # may not be square
    data = nozero_log(data, np.log2)
    vals = np.array([i for d in data for i in d])
    vals = vals[np.isfinite(vals)]

    mindata = np.nanmin(vals)
    maxdata = np.nanmax(vals)
    diff = maxdata - mindata
    posI = 0.01 if not clim else (float(clim[0]) / diff) if clim[0] != None else 0.01
    posF = 1.0  if not clim else (float(clim[1]) / diff) if clim[1] != None else 1.0
    if cmap == 'tadbit':
        cuts = perc
        cdict = {'red'  : [(0.0,  0.0, 0.0)],
                 'green': [(0.0,  0.0, 0.0)],
                 'blue' : [(0.0,  0.5, 0.5)]}
        prev_pos  = 0
        median = (np.median(vals) - mindata) / diff
        for prc in np.linspace(posI, median, cuts / 2, endpoint=False):
            try:
                pos = (np.percentile(vals, prc * 100.) - mindata) / diff
                prc = ((prc - posI) / (median - posI)) + 1. / cuts
            except ValueError:
                pos = prc = 0
            if prev_pos >= pos:
                continue
            cdict['red'  ].append([pos, prc, prc])
            cdict['green'].append([pos, prc, prc])
            cdict['blue' ].append([pos, 1, 1])
            prev_pos  = pos
        for prc in np.linspace(median + 1. / cuts, posF, cuts / 2, endpoint=False):
            try:
                pos = (np.percentile(vals, prc * 100.) - mindata) / diff
                prc = ((prc - median) / (posF - median))
            except ValueError:
                pos = prc = 0
            if prev_pos >= pos:
                continue
            cdict['red'  ].append([pos, 1.0, 1.0])
            cdict['green'].append([pos, 1 - prc, 1 - prc])
            cdict['blue' ].append([pos, 1 - prc, 1 - prc])
            prev_pos  = pos
        pos = (np.percentile(vals ,97.) - mindata) / diff
        cdict['red'  ].append([pos, 0.1, 0.1])
        cdict['green'].append([pos, 0, 0])
        cdict['blue' ].append([pos, 0, 0])

        cdict['red'  ].append([1.0, 1, 1])
        cdict['green'].append([1.0, 1, 1])
        cdict['blue' ].append([1.0, 0, 0])
        cmap  = LinearSegmentedColormap(cmap, cdict)
        clim = None
    else:
        cmap = plt.get_cmap(cmap)
    cmap.set_bad('darkgrey', 1)

    ax1.imshow(data, interpolation='none',
               cmap=cmap, vmin=clim[0] if clim else None, vmax=clim[1] if clim else None)
    size1 = len(data)
    size2 = len(data[0])
    if size1 == size2:
        for i in xrange(size1):
            for j in xrange(i, size2):
                if np.isnan(data[i][j]):
                    data[i][j] = 0
                    data[j][i] = 0
    else:
        for i in xrange(size1):
            for j in xrange(size2):
                if np.isnan(data[i][j]):
                    data[i][j] = 0
            #data[j][i] = data[i][j]
    try:
        evals, evect = eigh(data)
        sort_perm = evals.argsort()
        evect = evect[sort_perm]
    except:
        evals, evect = None, None
    data = [i for d in data for i in d if not np.isnan(i)]
    gradient = np.linspace(np.nanmin(data),
                           np.nanmax(data), max(size1, size2))
    gradient = np.vstack((gradient, gradient))
    h  = ax2.hist(data, color='darkgrey', linewidth=2,
                  bins=20, histtype='step', normed=True)
    _  = ax2.imshow(gradient, aspect='auto', cmap=cmap,
                    extent=(np.nanmin(data), np.nanmax(data) , 0, max(h[0])))
    if genome_seq:
        for crm in genome_seq:
            ax1.vlines([cumcs[crm][0]-.5, cumcs[crm][1]-.5], cumcs[crm][0]-.5, cumcs[crm][1]-.5,
                       color='w', linestyle='-', linewidth=1, alpha=1)
            ax1.hlines([cumcs[crm][1]-.5, cumcs[crm][0]-.5], cumcs[crm][0]-.5, cumcs[crm][1]-.5,
                       color='w', linestyle='-', linewidth=1, alpha=1)
            ax1.vlines([cumcs[crm][0]-.5, cumcs[crm][1]-.5], cumcs[crm][0]-.5, cumcs[crm][1]-.5,
                       color='k', linestyle='--')
            ax1.hlines([cumcs[crm][1]-.5, cumcs[crm][0]-.5], cumcs[crm][0]-.5, cumcs[crm][1]-.5,
                       color='k', linestyle='--')
        if not one:
            vals = [0]
            keys = ['']
            for crm in genome_seq:
                vals.append(cumcs[crm][0])
                keys.append(crm)
            vals.append(cumcs[crm][1])
            ax1.set_yticks(vals)
            ax1.set_yticklabels('')
            ax1.set_yticks([float(vals[i]+vals[i+1])/2
                            for i in xrange(len(vals) - 1)], minor=True)
            ax1.set_yticklabels(keys, minor=True)
            for t in ax1.yaxis.get_minor_ticks():
                t.tick1On = False
                t.tick2On = False
    # totaloridata = ''.join([j + ('' if (i+1)%3 else ',') for i, j in enumerate(str(totaloridata)[::-1])])[::-1].strip(',')
    # minoridata = ''.join([j + ('' if (i+1)%3 else ',') for i, j   in enumerate(str(minoridata)[::-1])])[::-1].strip(',')
    # maxoridata = ''.join([j + ('' if (i+1)%3 else ',') for i, j   in enumerate(str(maxoridata)[::-1])])[::-1].strip(',')
    plt.figtext(0.05,0.25, ''.join([
        (name + '\n') if name else '',
        'Number of interactions: %s\n' % str(totaloridata),
        ('' if np.isnan(cistrans) else
         ('Percentage of cis interactions: %.0f%%\n' % (cistrans*100))),
        'Min interactions: %s\n' % (minoridata),
        'Max interactions: %s\n' % (maxoridata)]))
    ax2.set_xlim((np.nanmin(data), np.nanmax(data)))
    ax2.set_ylim((0, max(h[0])))
    ax1.set_xlim ((-0.5, size1 - .5))
    ax1.set_ylim ((-0.5, size2 - .5))
    ax2.set_xlabel('log interaction count')
    # we reduce the number of dots displayed.... we just want to see the shape
    subdata = np.array(list(set([float(int(d*100))/100 for d in data])))
    try:
        normfit = sc_norm.pdf(subdata, np.nanmean(data), np.nanstd(data))
    except AttributeError:
        normfit = sc_norm.pdf(subdata, np.mean(data), np.std(data))
    ax2.plot(subdata, normfit, 'w.', markersize=2.5, alpha=.4)
    ax2.plot(subdata, normfit, 'k.', markersize=1.5, alpha=1)
    ax2.set_title('skew: %.3f, kurtosis: %.3f' % (skew(data),
                                                   kurtosis(data)))
    try: 
        ax4.vlines(range(size1), 0, evect[:,-1], color='k')
    except (TypeError, IndexError):
        pass
    ax4.hlines(0, 0, size2, color='red')
    ax4.set_ylabel('E1')
    ax4.set_yticklabels([])
    try:
        ax5.vlines(range(size1), 0, evect[:,-2], color='k')
    except (TypeError, IndexError):
        pass
    ax5.hlines(0, 0, size2, color='red')
    ax5.set_ylabel('E2')
    ax5.set_yticklabels([])
    try:
        ax6.vlines(range(size1), 0, evect[:,-3], color='k')
    except (TypeError, IndexError):
        pass
    ax6.hlines(0, 0, size2, color='red')
    ax6.set_ylabel('E3')
    ax6.set_yticklabels([])
    xticklabels = ax4.get_xticklabels() + ax5.get_xticklabels() + ax6.get_xticklabels()
    plt.setp(xticklabels, visible=False)
    if savefig:
        tadbit_savefig(savefig)
    elif show:
        plt.show()
    plt.close('all')


def plot_distance_vs_interactions(data, min_diff=1, max_diff=1000, show=False,
                                  genome_seq=None, resolution=None, axe=None,
                                  savefig=None, normalized=False,
                                  plot_each_cell=False):
    """
    :param data: input file name, or HiC_data object or list of lists
    :param 10 min_diff: lower limit (in number of bins)
    :param 1000 max_diff: upper limit (in number of bins) to look for
    :param 100 resolution: group reads that are closer than this resolution
       parameter
    :param_hash False plot_each_cell: if false, only the mean distances by bin
       will be represented, otherwise each pair of interactions will be plotted.
    :param None axe: a matplotlib.axes.Axes object to define the plot
       appearance
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).

    :returns: slope, intercept and R square of each of the 3 correlations
    """
    resolution = resolution or 1
    if isinstance(data, str):
        dist_intr = dict([(i, {})
                          for i in xrange(min_diff, max_diff)])
        fhandler = open(data)
        line = fhandler.next()
        while line.startswith('#'):
            line = fhandler.next()
        try:
            while True:
                _, cr1, ps1, _, _, _, _, cr2, ps2, _ = line.split('\t', 9)
                if cr1 != cr2:
                    line = fhandler.next()
                    continue
                diff = abs(int(ps1)  / resolution - int(ps2) / resolution)
                if max_diff > diff >= min_diff:
                    try:
                        dist_intr[diff][int(ps1) / resolution] += 1.
                    except KeyError:
                        dist_intr[diff][int(ps1) / resolution] = 1.
                line = fhandler.next()
        except StopIteration:
            pass
        fhandler.close()
        for diff in dist_intr:
            dist_intr[diff] = [dist_intr[diff].get(k, 0)
                               for k in xrange(max(dist_intr[diff]) - diff)]
    elif isinstance(data, HiC_data):
        dist_intr = dict([(i, []) for i in xrange(min_diff, max_diff)])
        if normalized:
            get_data = lambda x, y: data[x, y] / data.bias[x] / data.bias[y]
        else:
            get_data = lambda x, y: data[x, y]
        max_diff = min(len(data), max_diff)
        if data.section_pos:
            for crm in data.section_pos:
                for diff in xrange(min_diff, min(
                    (max_diff, 1 + data.chromosomes[crm]))):
                    for i in xrange(data.section_pos[crm][0],
                                    data.section_pos[crm][1] - diff):
                        dist_intr[diff].append(get_data(i, i + diff))
        else:
            for diff in xrange(min_diff, max_diff):
                for i in xrange(len(data) - diff):
                    if not np.isnan(data[i, i + diff]):
                        dist_intr[diff].append(get_data(i, diff))
    else:
        dist_intr = dict([(i, []) for i in xrange(min_diff, max_diff)])
        if genome_seq:
            max_diff = min(max(genome_seq.values()), max_diff)
            cnt = 0
            for crm in genome_seq:
                for diff in xrange(min_diff, min(
                    (max_diff, genome_seq[crm]))):
                    for i in xrange(cnt, cnt + genome_seq[crm] - diff):
                        if not np.isnan(data[i][i + diff]):
                            dist_intr[diff].append(data[i][i + diff])
                cnt += genome_seq[crm]
        else:
            max_diff = min(len(data), max_diff)
            for diff in xrange(min_diff, max_diff):
                for i in xrange(len(data) - diff):
                    if not np.isnan(data[i][i + diff]):
                        dist_intr[diff].append(data[i][i + diff])
    if not axe:
        fig=plt.figure()
        axe = fig.add_subplot(111)
    # remove last part of the plot in case no interaction is count... reduce max_dist
    for diff in xrange(max_diff - 1, min_diff, -1):
        try:
            if not dist_intr[diff]:
                del(dist_intr[diff])
                max_diff -=1
                continue
        except KeyError:
            max_diff -=1
            continue
        break
    # get_cmap the mean values perc bins
    mean_intr = dict([(i, float(sum(dist_intr[i])) / len(dist_intr[i]))
                      for i in dist_intr if len(dist_intr[i])])
    if plot_each_cell:
        xp, yp = [], []
        for x, y in sorted(dist_intr.items(), key=lambda x:x[0]):
            xp.extend([x] * len(y))
            yp.extend(y)
        x = []
        y = []
        for k in xrange(len(xp)):
            if yp[k]:
                x.append(xp[k])
                y.append(yp[k])
        axe.plot(x, y, color='grey', marker='.', alpha=0.1, ms=1,
                 linestyle='None')
    xp, yp = zip(*sorted(mean_intr.items(), key=lambda x:x[0]))
    x = []
    y = []
    for k in xrange(len(xp)):
        if yp[k]:
            x.append(xp[k])
            y.append(yp[k])
    axe.plot(x, y, 'k.')
    best = (float('-inf'), 0, 0, 0, 0, 0, 0, 0, 0, 0)
    logx = np.log(x)
    logy = np.log(y)
    ntries = 100
    # set k for better fit
    # for k in xrange(1, ntries/5, ntries/5/5):
    if resolution == 1:
        k = 1
        for i in xrange(3, ntries-2-k):
            v1 = i * len(x) / ntries
            try:
                a1, b1, r21, _, _ = linregress(logx[ :v1], logy[ :v1])
            except ValueError:
                a1 = b1 = r21 = 0
            r21 *= r21
            for j in xrange(i + 1 + k, ntries - 2 - k):
                v2 = j * len(x) / ntries
                try:
                    a2, b2, r22, _, _ = linregress(logx[v1+k:v2], logy[v1+k:v2])
                    a3, b3, r23, _, _ = linregress(logx[v2+k:  ], logy[v2+k: ])
                except ValueError:
                    a2 = b2 = r22 = 0
                    a3 = b3 = r23 = 0
                r2 = r21 + r22**2 + r23**2
                if r2 > best[0]:
                    best = (r2, v1, v2, a1, a2, a3,
                            b1, b2, b3, k)
        # plot line of best fit
        (v1, v2, 
         a1, a2, a3,
         b1, b2, b3, k) = best[1:]
        yfit1 = lambda xx: np.exp(b1 + a1*np.array (np.log(xx)))
        yfit2 = lambda xx: np.exp(b2 + a2*np.array (np.log(xx)))
        yfit3 = lambda xx: np.exp(b3 + a3*np.array (np.log(xx)))
        axe.plot(x[  :v1], yfit1(x[  :v1] ), color= 'yellow', lw=2,
                 label = r'$\alpha_{%s}=%.2f$' % (
                     '0-0.7 \mathrm{ Mb}' if resolution != 1 else '1', a1))
                 #label = r'$\alpha_1=%.2f$ (0-%d)' % (a1, x[v1]))
        axe.plot(x[v1+k:v2], yfit2(x[v1+k:v2]),  color= 'orange', lw=2,
                 label = r'$\alpha_{%s}=%.2f$' % (
                     '0.7-10 \mathrm{ Mb}' if resolution != 1 else '2', a2))
                 # label = r'$\alpha_2=%.2f$ (%d-%d)' % (a2, x[v1], x[v2]))
        axe.plot(x[v2+k:  ], yfit3(x[v2+k:  ] ), color= 'red'   , lw=2,
                 label = r'$\alpha_{%s}=%.2f$' % (
                     '10 \mathrm{ Mb}-\infty' if resolution != 1 else '3', a3))
                 # label = r'$\alpha_3=%.2f$ (%d-$\infty$)' % (a3, x[v2+k]))
    else:
        # from 0.7 Mb
        v1 = 700000   / resolution
        # to 10 Mb
        v2 = 10000000 / resolution
        try:
            a1, b1, r21, _, _ = linregress(logx[  :v1], logy[  :v1])
        except ValueError:
            a1, b1, r21 = 0, 0, 0
        try:
            a2, b2, r22, _, _ = linregress(logx[v1:v2], logy[v1:v2])
        except ValueError:
            a2, b2, r22 = 0, 0, 0
        try:
            a3, b3, r23, _, _ = linregress(logx[v2:  ], logy[v2:  ])
        except ValueError:
            a3, b3, r23 = 0, 0, 0
        yfit1 = lambda xx: np.exp(b1 + a1*np.array (np.log(xx)))
        yfit2 = lambda xx: np.exp(b2 + a2*np.array (np.log(xx)))
        yfit3 = lambda xx: np.exp(b3 + a3*np.array (np.log(xx)))
        axe.plot(x[  :v1], yfit1(x[  :v1] ), color= 'yellow', lw=2,
                 label = r'$\alpha_{%s}=%.2f$' % (
                     '0-0.7 \mathrm{ Mb}' if resolution != 1 else '1', a1))
                 #label = r'$\alpha_1=%.2f$ (0-%d)' % (a1, x[v1]))
        axe.plot(x[v1:v2], yfit2(x[v1:v2]),  color= 'orange', lw=2,
                 label = r'$\alpha_{%s}=%.2f$' % (
                     '0.7-10 \mathrm{ Mb}' if resolution != 1 else '2', a2))
                 # label = r'$\alpha_2=%.2f$ (%d-%d)' % (a2, x[v1], x[v2]))
        axe.plot(x[v2:  ], yfit3(x[v2:  ] ), color= 'red'   , lw=2,
                 label = r'$\alpha_{%s}=%.2f$' % (
                     '10 \mathrm{ Mb}-\infty' if resolution != 1 else '3', a3))
                 # label = r'$\alpha_3=%.2f$ (%d-$\infty$)' % (a3, x[v2+k]))
    axe.set_ylabel('Log interaction count')
    axe.set_xlabel('Log genomic distance (resolution: %s)' % nicer(resolution))
    axe.legend(loc='lower left', frameon=False)
    axe.set_xscale('log')
    axe.set_yscale('log')
    axe.set_xlim((min_diff, max_diff))
    try:
        axe.set_ylim((0, max(y)))
    except ValueError:
        pass
    if savefig:
        tadbit_savefig(savefig)
        plt.close('all')
    elif show==True:
        plt.show()
        plt.close('all')
    return (a1, b1, r21), (a2, b2, r22), (a3, b3, r23)

def plot_iterative_mapping(fnam1, fnam2, total_reads=None, axe=None, savefig=None):
    """
    :param fnam: input file name
    :param total_reads: total number of reads in the initial FASTQ file
    :param None axe: a matplotlib.axes.Axes object to define the plot
       appearance
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :returns: a dictionary with the number of reads per mapped length
    """
    count_by_len = {}
    total_reads = total_reads or 1
    if not axe:
        fig=plt.figure()
        _ = fig.add_subplot(111)
    colors = ['olive', 'darkcyan']
    iteration = False
    for i, fnam in enumerate([fnam1, fnam2]):
        fhandler = open(fnam)
        line = fhandler.next()
        count_by_len[i] = {}
        while line.startswith('#'):
            if line.startswith('# MAPPED '):
                itr, num = line.split()[2:]
                count_by_len[i][int(itr)] = int(num)
            line = fhandler.next()
        if not count_by_len[i]:
            iteration = True
            try:
                while True:
                    _, length, _, _ = line.rsplit('\t', 3)
                    try:
                        count_by_len[i][int(length)] += 1
                    except KeyError:
                        count_by_len[i][int(length)] = 1
                    line = fhandler.next()
            except StopIteration:
                pass
        fhandler.close()
        lengths = sorted(count_by_len[i].keys())
        for k in lengths[::-1]:
            count_by_len[i][k] += sum([count_by_len[i][j]
                                       for j in lengths if j < k])
        plt.plot(lengths, [float(count_by_len[i][l]) / total_reads
                           for l in lengths],
                 label='read' + str(i + 1), linewidth=2, color=colors[i])
    if iteration:
        plt.xlabel('read length (bp)')
    else:
        plt.xlabel('Iteration number')
    if total_reads != 1:
        plt.ylabel('Proportion of mapped reads')
    else:
        plt.ylabel('Number of mapped reads')
    plt.legend(loc=4)
    if savefig:
        tadbit_savefig(savefig)
    elif not axe:
        plt.show()
    plt.close('all')
    return count_by_len


def insert_sizes(fnam, savefig=None, nreads=None, max_size=99.9, axe=None,
                 show=False, xlog=False, stats=('median', 'perc_max')):
    """
    Plots the distribution of dangling-ends lengths
    :param fnam: input file name
    :param None savefig: path where to store the output images.
    :param 99.9 max_size: top percentage of distances to consider, within the
       top 0.01% are usually found very long outliers.
    :param False xlog: represent x axis in logarithmic scale
    :param ('median', 'perc_max') stats: returns this set of values calculated from the
       distribution of insert/fragment sizes. Possible values are:
        - 'median' median of the distribution
        - 'perc_max' percentil defined by the other parameter 'max_size'
        - 'first_deacay' starting from the median of the distribution to the
            first window where 10 consecutive insert sizes are counted less than
            a given value (this given value is equal to the sum of all
            sizes divided by 100 000)
        - 'MAD' Double Median Adjusted Deviation

    :returns: the median value and the percentile inputed as max_size.
    """
    distr = {}
    genome_seq = OrderedDict()
    pos = 0
    fhandler = open(fnam)
    for line in fhandler: 
        if line.startswith('#'):
            if line.startswith('# CRM '):
                crm, clen = line[6:].split('\t')
                genome_seq[crm] = int(clen)
        else:
            break
        pos += len(line)
    fhandler.seek(pos)
    des = []
    if nreads:
        nreads /= 2
    for line in fhandler:
        (crm1, pos1, dir1, _, re1, _,
         crm2, pos2, dir2, _, re2) = line.strip().split('\t')[1:12]
        if re1==re2 and crm1 == crm2 and dir1 != dir2:
            pos1, pos2 = int(pos1), int(pos2)
            if (pos2 > pos1) == int(dir1):
                des.append(abs(pos2 - pos1))
            if len(des) == nreads:
                break
    fhandler.close()
    max_perc = np.percentile(des, max_size)
    perc99   = np.percentile(des, 99)
    perc01   = np.percentile(des, 1)
    perc50   = np.percentile(des, 50)
    perc95   = np.percentile(des, 95)
    perc05   = np.percentile(des, 5)
    to_return = {'median': perc50}
    cutoff = len(des) / 100000.
    count  = 0
    for v in xrange(int(perc50), int(max(des))):
        if des.count(v) < cutoff:
            count += 1
        else:
            count = 0
        if count >= 10:
            to_return['first_decay'] = v - 10
            break
    else:
        raise Exception('ERROR: not found')
    to_return['perc_max'] = max_perc
    to_return['MAD'] = mad(des)
    if not savefig and not axe and not show:
        return [to_return[k] for k in stats]
    
    ax = setup_plot(axe, figsize=(10, 5.5))
    desapan = ax.axvspan(perc95, perc99, facecolor='darkolivegreen', alpha=.3,
                         label='1-99%% DEs\n(%.0f-%.0f nts)' % (perc01, perc99))
    ax.axvspan(perc01, perc05, facecolor='darkolivegreen', alpha=.3)
    desapan = ax.axvspan(perc05, perc95, facecolor='darkseagreen', alpha=.3,
                         label='5-95%% DEs\n(%.0f-%.0f nts)' % (perc05, perc95))
    deshist = ax.hist(des, bins=100, range=(0, max_perc),
                      alpha=.7, color='darkred', label='Dangling-ends')
    ylims   = ax.get_ylim()
    plots   = []
    ax.set_xlabel('Genomic distance between reads')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of dangling-ends ' +
                 'lenghts\n(median: %s, top %.1f%%, up to %0.f nts)' % (
                     perc50, max_size, max_perc))
    if xlog:
        ax.set_xscale('log')
    ax.set_xlim((50, max_perc))
    plt.subplots_adjust(left=0.1, right=0.75)
    ax.legend(bbox_to_anchor=(1.4, 1), frameon=False)
    if savefig:
        tadbit_savefig(savefig)
    elif show and not axe:
        plt.show()
    plt.close('all')
    return [to_return[k] for k in stats]


def plot_genomic_distribution(fnam, first_read=True, resolution=10000,
                              ylim=None, yscale=None, savefig=None, show=False,
                              savedata=None, chr_names=None, nreads=None):
    """
    :param fnam: input file name
    :param True first_read: uses first read.
    :param 100 resolution: group reads that are closer than this resolution
       parameter
    :param None ylim: a tuple of lower and upper bound for the y axis
    :param None yscale: if set_bad to "log" values will be represented in log2 
       scale
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param None savedata: path where to store the output read counts per bin.
    :param None chr_names: can pass a list of chromosome names in case only some
       them the need to be plotted (this option may last even more than default)
    
    """
    distr = {}
    idx1, idx2 = (1, 3) if first_read else (7, 9)
    genome_seq = OrderedDict()
    if chr_names:
        chr_names = set(chr_names)
        cond1 = lambda x: x not in chr_names
    else:
        cond1 = lambda x: False
    if nreads:
        cond2 = lambda x: x >= nreads
    else:
        cond2 = lambda x: False
    cond = lambda x, y: cond1(x) and cond2(y)
    count = 0
    pos = 0
    fhandler = open(fnam)
    for line in fhandler: 
        if line.startswith('#'):
            if line.startswith('# CRM '):
                crm, clen = line[6:].split('\t')
                genome_seq[crm] = int(clen)
        else:
            break
        pos += len(line)
    fhandler.seek(pos)
    for line in fhandler:
        crm, pos = line.strip().split('\t')[idx1:idx2]
        count += 1
        if cond(crm, count):
            line = fhandler.next()
            if cond2(count):
                break
            continue
        pos = int(pos) / resolution
        try:
            distr[crm][pos] += 1
        except KeyError:
            try:
                distr[crm][pos] = 1
            except KeyError:
                distr[crm] = {pos: 1}
    fhandler.close()
    if savefig or show:
        _ = plt.figure(figsize=(15, 1 + 3 * len(
            chr_names if chr_names else distr.keys())))

    max_y = max([max(distr[c].values()) for c in distr])
    max_x = max([len(distr[c].values()) for c in distr])
    ncrms = len(chr_names if chr_names else genome_seq if genome_seq else distr)
    data = {}
    for i, crm in enumerate(chr_names if chr_names else genome_seq
                            if genome_seq else distr):
        try:
            # data[crm] = [distr[crm].get(j, 0) for j in xrange(max(distr[crm]))]  # genome_seq[crm]
            data[crm] = [distr[crm].get(j, 0)
                         for j in xrange(genome_seq[crm] / resolution + 1)]
            if savefig or show:
                plt.subplot(ncrms, 1, i + 1)
                plt.plot(range(max(distr[crm])), data[crm],
                         color='red', lw=1.5, alpha=0.7)
                if yscale:
                    plt.yscale(yscale)
        except KeyError:
            pass
        if savefig or show:
            if ylim:
                plt.vlines(genome_seq[crm] / resolution, ylim[0], ylim[1])
            else:
                plt.vlines(genome_seq[crm] / resolution, 0, max_y)
            plt.xlim((0, max_x))
            plt.ylim(ylim or (0, max_y))
            plt.title(crm)

    if savefig:
        tadbit_savefig(savefig)
        if not show:
            plt.close('all')
    elif show:
        plt.show()

    if savedata:
        out = open(savedata, 'w')
        out.write('# CRM\tstart-end\tcount\n')
        out.write('\n'.join('%s\t%d-%d\t%d' % (c, (i * resolution) + 1,
                                               ((i + 1) * resolution), v)
                            for c in data for i, v in enumerate(data[c])))
        out.write('\n')
        out.close()

def correlate_matrices(hic_data1, hic_data2, max_dist=10, intra=False, axe=None,
                       savefig=None, show=False, savedata=None,
                       normalized=False, remove_bad_columns=True, **kwargs):
    """
    Compare the iteractions of two Hi-C matrices at a given distance,
    with spearman rank correlation

    :param hic_data1: Hi-C-data object
    :param hic_data2: Hi-C-data object
    :param 1 resolution: to be used for scaling the plot
    :param 10 max_dist: maximum distance from diagonal (e.g. 10 mean we will
       not look further than 10 times the resolution)
    :param None savefig: path to save the plot
    :param False intra: only takes into account intra-chromosomal contacts
    :param False show: displays the plot
    :param False normalized: use normalized data
    :param True remove_bads: computes the union of bad columns between samples
       and exclude them from the comparison

    :returns: list of correlations and list of genomic distances
    """
    corrs = []
    dists = []

    if normalized:
        get_the_guy1 = lambda i, j: (hic_data1[j, i] / hic_data1.bias[i] /
                                     hic_data1.bias[j])
        get_the_guy2 = lambda i, j: (hic_data2[j, i] / hic_data2.bias[i] /
                                     hic_data2.bias[j])
    else:
        get_the_guy1 = lambda i, j: hic_data1[j, i]
        get_the_guy2 = lambda i, j: hic_data2[j, i]
    
    if remove_bad_columns:
        # union of bad columns
        bads = hic_data1.bads.copy()
        bads.update(hic_data2.bads)

    if (intra and hic_data1.sections and hic_data2.sections and 
        hic_data1.sections == hic_data2.sections):
        for dist in xrange(1, max_dist + 1):
            diag1 = []
            diag2 = []
            for crm in hic_data1.section_pos:
                for j in xrange(hic_data1.section_pos[crm][0],
                                hic_data1.section_pos[crm][1] - dist):
                    i = j + dist
                    if j in bads or i in bads:
                        continue
                    diag1.append(get_the_guy1(i, j))
                    diag2.append(get_the_guy2(i, j))
            corrs.append(spearmanr(diag1, diag2)[0])
            dists.append(dist)
    else:
        if intra:
            warn('WARNING: hic_dta does not contain chromosome coordinates, ' +
                 'intra set to False')
        for dist in xrange(1, max_dist + 1):
            diag1 = []
            diag2 = []
            for j in xrange(len(hic_data1) - dist):
                i = j + dist
                if j in bads or i in bads:
                    continue
                diag1.append(get_the_guy1(i, j))
                diag2.append(get_the_guy2(i, j))
            corrs.append(spearmanr(diag1, diag2)[0])
            dists.append(dist)
    if show or savefig or axe:
        if not axe:
            fig = plt.figure()
            axe = fig.add_subplot(111)
            given_axe = False
        else:
            given_axe = True
        axe.plot(dists, corrs, color='orange', linewidth=3, alpha=.8)
        axe.set_xlabel('Genomic distance in bins')
        axe.set_ylabel('Spearman rank correlation')
        axe.set_xlim((0, dists[-1]))
        if savefig:
            tadbit_savefig(savefig)
        if show:
            plt.show()
        if not given_axe:
            plt.close('all')
    if savedata:
        out = open(savedata, 'w')
        out.write('# genomic distance\tSpearman rank correlation\n')
        for i in xrange(len(corrs)):
            out.write('%s\t%s\n' % (dists[i], corrs[i]))
        out.close()
    if kwargs.get('get_bads', False):
        return corrs, dists, bads
    else:
        return corrs, dists

def eig_correlate_matrices(hic_data1, hic_data2, nvect=6, normalized=False, 
                           savefig=None, show=False, savedata=None,
                           remove_bad_columns=True, **kwargs):
    """
    Compare the iteractions of two Hi-C matrices using their 6 first
    eigenvectors, with pearson correlation

    :param hic_data1: Hi-C-data object
    :param hic_data2: Hi-C-data object
    :param 6 nvect: number of eigenvectors to compare
    :param None savefig: path to save the plot
    :param False show: displays the plot
    :param False normalized: use normalized data
    :param True remove_bads: computes the union of bad columns between samples
       and exclude them from the comparison
    :param kwargs: any argument to pass to matplotlib imshow function

    :returns: matrix of correlations
    """
    data1 = hic_data1.get_matrix(normalized=normalized)
    data2 = hic_data2.get_matrix(normalized=normalized)
    ## reduce matrices to remove bad columns
    if remove_bad_columns:
        # union of bad columns
        bads = hic_data1.bads.copy()
        bads.update(hic_data2.bads)
        # remove them form both matrices
        for bad in sorted(bads, reverse=True):
            del(data1[bad])
            del(data2[bad])
            for i in xrange(len(data1)):
                _ = data1[i].pop(bad)
                _ = data2[i].pop(bad)
    # get the log
    size = len(data1)
    data1 = nozero_log(data1, np.log2)
    data2 = nozero_log(data2, np.log2)
    # get the eigenvectors
    ev1, evect1 = eigh(data1)
    ev2, evect2 = eigh(data2)
    corr = [[0 for _ in xrange(nvect)] for _ in xrange(nvect)]
    # sort eigenvectors according to their eigenvalues => first is last!!
    sort_perm = ev1.argsort()
    ev1.sort()
    evect1 = evect1[sort_perm]
    sort_perm = ev2.argsort()
    ev2.sort()
    evect2 = evect2[sort_perm]
    # calculate Pearson correlation
    for i in xrange(nvect):
        for j in xrange(nvect):
            corr[i][j] = abs(pearsonr(evect1[:,-i-1],
                                      evect2[:,-j-1])[0])
    # plot
    axe    = plt.axes([0.1, 0.1, 0.6, 0.8])
    cbaxes = plt.axes([0.85, 0.1, 0.03, 0.8])
    if show or savefig:
        im = axe.imshow(corr, interpolation="nearest",origin='lower', **kwargs)
        axe.set_xlabel('Eigen Vectors exp. 1')
        axe.set_ylabel('Eigen Vectors exp. 2')
        axe.set_xticks(range(nvect))
        axe.set_yticks(range(nvect))
        axe.set_xticklabels(range(1, nvect + 2))
        axe.set_yticklabels(range(1, nvect + 2))
        axe.xaxis.set_tick_params(length=0, width=0)
        axe.yaxis.set_tick_params(length=0, width=0)
        
        cbar = plt.colorbar(im, cax = cbaxes )
        cbar.ax.set_ylabel('Pearson correlation', rotation=90*3,
                           verticalalignment='bottom')
        axe2 = axe.twinx()
        axe2.set_yticks(range(nvect))
        axe2.set_yticklabels(['%.1f' % (e) for e in ev2[-nvect:][::-1]])
        axe2.set_ylabel('corresponding Eigen Values exp. 2', rotation=90*3,
                        verticalalignment='bottom')
        axe2.set_ylim((-0.5, nvect - 0.5))
        axe2.yaxis.set_tick_params(length=0, width=0)
        
        axe3 = axe.twiny()
        axe3.set_xticks(range(nvect))
        axe3.set_xticklabels(['%.1f' % (e) for e in ev1[-nvect:][::-1]])
        axe3.set_xlabel('corresponding Eigen Values exp. 1')
        axe3.set_xlim((-0.5, nvect - 0.5))
        axe3.xaxis.set_tick_params(length=0, width=0)
        
        axe.set_ylim((-0.5, nvect - 0.5))
        axe.set_xlim((-0.5, nvect - 0.5))
        if savefig:
            tadbit_savefig(savefig)
        if show:
            plt.show()
        plt.close('all')

    if savedata:
        out = open(savedata, 'w')
        out.write('# ' + '\t'.join(['Eigen Vector %s'% i
                                    for i in xrange(nvect)]) + '\n')
        for i in xrange(nvect):
            out.write('\t'.join([str(corr[i][j])
                                 for j in xrange(nvect)]) + '\n')
        out.close()
    if kwargs.get('get_bads', False):
        return corr, bads
    else:
        return corr

def plot_rsite_reads_distribution(reads_file, outprefix, window=20,
        maxdist=1000):
    de_right={}
    de_left={}
    print "process reads"
    fl=open(reads_file)
    while True:
        line=fl.next()
        if not line.startswith('#'):
            break
    nreads=0
    try:
        while True:
            nreads += 1
            if nreads % 1000000 == 0:
                print nreads
            try:
                _, n1, sb1, sd1, l1, ru1, rd1, n2, sb2, sd2, l2, ru2, rd2\
                        = line.split()
                sb1, sd1, l1, ru1, rd1, sb2, sd2, l2, ru2, rd2 = \
                        map(int, [sb1, sd1, l1, ru1, rd1, sb2, sd2, l2,
                            ru2, rd2])
            except ValueError:
                print line
                raise ValueError("line is not the right format!")
            if n1 != n2:
                line=fl.next()
                continue
            #read1 ahead of read2
            if sb1 > sb2:
                sb1, sd1, l1, ru1, rd1, sb2, sd2, l2, ru2, rd2 = \
                    sb2, sd2, l2, ru2, rd2, sb1, sd1, l1, ru1, rd1
            #direction always -> <-
            if not (sd1 == 1 and sd2 == 0):
                line=fl.next()
                continue
            #close to the diagonal
            if sb2-sb1 > maxdist:
                line=fl.next()
                continue
            #close to RE 1
            if abs(sb1-ru1) < abs(sb1-rd1):
                rc1=ru1
            else:
                rc1=rd1
            pos=sb1-rc1
            if abs(pos)<=window:
                if not pos in de_right:
                    de_right[pos]=0
                de_right[pos]+=1
            #close to RE 2
            if abs(sb2-ru2) < abs(sb2-rd2):
                rc2=ru2
            else:
                rc2=rd2
            pos=sb2-rc2
            if abs(pos)<=window:
                if not pos in de_left:
                    de_left[pos]=0
                de_left[pos]+=1
            line=fl.next()
    except StopIteration:
        pass
    print "   finished processing {} reads".format(nreads)

    #transform to arrays
    ind = range(-window,window+1)
    de_r = map(lambda x:de_right.get(x,0), ind)
    de_l = map(lambda x:de_left.get(x,0), ind)

    #write to files
    print "write to files"
    fl=open(outprefix+'_count.dat','w')
    fl.write('#dist\tX~~\t~~X\n')
    for i,j,k in zip(ind,de_r, de_l):
        fl.write('{}\t{}\t{}\n'.format(i, j, k))

    #write plot
    rcParams.update({'font.size': 10})
    pp = PdfPages(outprefix+'_plot.pdf')
    ind = np.array(ind)
    width = 1
    pr = plt.bar(ind-0.5, de_r, width, color='r')
    pl = plt.bar(ind-0.5, de_l, width, bottom=de_r, color='b')
    plt.ylabel("Count")
    plt.title("Histogram of counts around cut site")
    plt.xticks(ind[::2], rotation="vertical")
    plt.legend((pl[0], pr[0]), ("~~X", "X~~")) 
    plt.gca().set_xlim([-window-1,window+1])
    pp.savefig()
    pp.close()


def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def plot_diagonal_distributions(reads_file, outprefix, ma_window=20,
        maxdist=800, de_left=[-2,3], de_right=[0,5]):
    rbreaks={}
    rejoined={}
    des={}
    print "process reads"
    fl=open(reads_file)
    while True:
        line=fl.next()
        if not line.startswith('#'):
            break
    nreads=0
    try:
        while True:
            nreads += 1
            if nreads % 1000000 == 0:
                print nreads
            try:
                _, n1, sb1, sd1, _, ru1, rd1, n2, sb2, sd2, _, ru2, rd2\
                        = line.split()
                sb1, sd1, ru1, rd1, sb2, sd2, ru2, rd2 = \
                        map(int, [sb1, sd1, ru1, rd1, sb2, sd2, ru2, rd2])
            except ValueError:
                print line
                raise ValueError("line is not the right format!")
            if n1 != n2:
                line=fl.next()
                continue
            #read1 ahead of read2
            if sb1 > sb2:
                sb1, sd1, ru1, rd1, sb2, sd2, ru2, rd2 = \
                    sb2, sd2, ru2, rd2, sb1, sd1, ru1, rd1
            #direction always -> <-
            if not (sd1 == 1 and sd2 == 0):
                line=fl.next()
                continue
            mollen = sb2-sb1
            if mollen > maxdist:
                line=fl.next()
                continue
            #DE1
            if abs(sb1-ru1) < abs(sb1-rd1):
                rc1=ru1
            else:
                rc1=rd1
            pos=sb1-rc1
            if pos in de_right:
                if not mollen in des:
                    des[mollen]=0
                des[mollen]+=1
                line=fl.next()
                continue
            #DE2
            if abs(sb2-ru2) < abs(sb2-rd2):
                rc2=ru2
            else:
                rc2=rd2
            pos=sb2-rc2
            if pos in de_left:
                if not mollen in des:
                    des[mollen]=0
                des[mollen]+=1
                line=fl.next()
                continue
            #random: map on same fragment
            if rd1 == rd2:
                if not mollen in rbreaks:
                    rbreaks[mollen]=0
                rbreaks[mollen]+=1
                line=fl.next()
                continue
            #rejoined ends
            if not mollen in rejoined:
                rejoined[mollen]=0
            rejoined[mollen]+=1
            line=fl.next()
    except StopIteration:
        pass
    print "   finished processing {} reads".format(nreads)

    #transform to arrays
    maxlen = max(max(rejoined),max(des),max(rbreaks))
    ind = range(1,maxlen+1)
    des = map(lambda x:des.get(x,0), ind)
    rbreaks = map(lambda x:rbreaks.get(x,0), ind)
    rejoined = map(lambda x:rejoined.get(x,0), ind)
    #reweight corner for rejoined
    rejoined = map(lambda x: x**.5 * rejoined[x-1]/x, ind)

    #write to files
    print "write to files"
    fl=open(outprefix+'_count.dat','w')
    fl.write('#dist\trbreaks\tdes\trejoined\n')
    for i,j,k,l in zip(ind,rbreaks,des,rejoined):
        fl.write('{}\t{}\t{}\t{}\n'.format(i, j, k, l))

    #transform data a bit more
    ind, des, rbreaks, rejoined = \
            map(lambda x: moving_average(np.array(x), ma_window),
                    [ind, des, rbreaks, rejoined])
    des, rbreaks, rejoined = map(lambda x:x/float(x.sum()),
            [des, rbreaks, rejoined])
    np.insert(ind,0,0)
    np.insert(des,0,0)
    np.insert(rbreaks,0,0)
    np.insert(rejoined,0,0)

    #write plot
    pp = PdfPages(outprefix+'_plot.pdf')
    rcParams.update({'font.size': 10})
    pde = plt.fill_between(ind, des, 0, color='r', alpha=0.5)
    prb = plt.fill_between(ind, rbreaks, 0, color='b', alpha=0.5)
    prj = plt.fill_between(ind, rejoined, 0, color='y', alpha=0.5)
    plt.ylabel("Normalized count")
    plt.ylabel("Putative DNA molecule length")
    plt.title("Histogram of counts close to the diagonal")
    #plt.xticks(ind[::10], rotation="vertical")
    plt.legend((prb, pde, prj), ("Random breaks", "Dangling ends",
        "Rejoined"))
    plt.gca().set_xlim([0,maxlen])
    pp.savefig()
    pp.close()


