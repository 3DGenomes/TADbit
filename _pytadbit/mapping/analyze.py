"""
18 Nov 2014


"""
from pytadbit.utils.extraviews   import tadbit_savefig, setup_plot
from pytadbit.utils.tadmaths     import nozero_log_matrix as nozero_log
from pytadbit.parsers.hic_parser import HiC_data
from warnings                    import warn
from collections                 import OrderedDict
from pytadbit.parsers.hic_parser import load_hic_data_from_reads
from pytadbit.utils.extraviews   import nicer
from scipy.stats                 import norm as sc_norm, skew, kurtosis
from scipy.stats                 import pearsonr, spearmanr, linregress
from numpy.linalg                import eigh
import os
import numpy as np

try:
    from matplotlib import pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:
    warn('matplotlib not found\n')


def hic_map(data, resolution=None, normalized=False, masked=None,
            by_chrom=False, savefig=None, show=False, savedata=None,
            focus=None, clim=None, cmap='tadbit', pdf=False, decay=True,
            perc=10, name=None, **kwargs):
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
    :param None clim: cutoff for the upper and lower bound in the coloring scale
       of the heatmap
    :param False pdf: when using the bny_chrom option, to specify the format of
       the stored images
    :param Reds cmap: color map to be used for the heatmap
    """
    if isinstance(data, str):
        data = load_hic_data_from_reads(data, resolution=resolution, **kwargs)
        if not kwargs.get('get_sections', True) and decay:
            warn('WARNING: not decay not available when get_sections is off.')
            decay = False
    hic_data = data
    resolution = data.resolution
    if hic_data.bads and not masked:
        masked = hic_data.bads
    # save and draw the data
    if by_chrom:
        if focus:
            raise Exception('Incompatible options focus and by_chrom\n')
        os.system('mkdir -p ' + (savedata if savedata else savefig))
        for i, crm1 in enumerate(hic_data.chromosomes):
            for crm2 in hic_data.chromosomes.keys()[i:]:
                if by_chrom == 'intra' and crm1 != crm2:
                    continue
                if by_chrom == 'inter' and crm1 == crm2:
                    continue
                subdata = hic_data.get_matrix(focus=(crm1, crm2), normalized=normalized)
                if savedata:
                    hic_data.write_matrix('%s/%s.mat' % (
                        savedata, '_'.join(set((crm1, crm2)))),
                                          focus=(crm1, crm2),
                                          normalized=normalized)
                if show or savefig:
                    draw_map(subdata, 
                             OrderedDict([(k, hic_data.chromosomes[k])
                                          for k in hic_data.chromosomes.keys()
                                          if k in [crm1, crm2]]),
                             hic_data.section_pos,
                             '%s/%s.%s' % (savefig,
                                           '_'.join(set((crm1, crm2))),
                                           'pdf' if pdf else 'png'),
                             show, one=True, clim=clim, cmap=cmap,
                             resolution=resolution, perc=perc,
                             name=name, cistrans=float('NaN'))
    else:
        if savedata:
            hic_data.write_matrix(savedata, focus=focus,
                                  normalized=normalized)
        if show or savefig:
            subdata = hic_data.get_matrix(focus=focus, normalized=normalized)
            if focus and masked:
                # rescale masked
                masked = dict([(m - focus[0], masked[m]) for m in masked])
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
                     clim=clim, cmap=cmap, resolution=resolution, perc=perc,
                     name=name, cistrans=float('NaN') if focus else
                     hic_data.cis_trans_ratio(kwargs.get('normalized', False),
                                              kwargs.get('exclude', None),
                                              kwargs.get('diagonal', False),
                                              kwargs.get('equals', None),
                                              kwargs.get('verbose', False)))


def draw_map(data, genome_seq, cumcs, savefig, show, resolution=None, one=False,
             clim=None, cmap='tadbit', decay=False, perc=10, name=None, cistrans=None):
    _ = plt.figure(figsize=(15.,12.5))
    ax1 = plt.axes([0.34, 0.08, 0.6, 0.7205])
    ax2 = plt.axes([0.07, 0.65, 0.21, 0.15])
    if decay:
        ax3 = plt.axes([0.07, 0.42, 0.21, 0.15])
        plot_distance_vs_interactions(data, genome_seq=genome_seq, axe=ax3,
                                      resolution=resolution)
    ax4 = plt.axes([0.34, 0.805, 0.6, 0.04], sharex=ax1)
    ax5 = plt.axes([0.34, 0.845, 0.6, 0.04], sharex=ax1)
    ax6 = plt.axes([0.34, 0.885, 0.6, 0.04], sharex=ax1)
    minoridata   = np.min(data)
    maxoridata   = np.max(data)
    totaloridata = np.nansum([data[i][j] for i in xrange(len(data))
                              for j in xrange(i, len(data))])
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
    size = len(data)
    for i in xrange(size):
        for j in xrange(i, size):
            if np.isnan(data[i][j]):
                data[i][j] = 0
                data[j][i] = 0
            #data[j][i] = data[i][j]
    evals, evect = eigh(data)
    sort_perm = evals.argsort()
    evect = evect[sort_perm]
    data = [i for d in data for i in d if not np.isnan(i)]
    gradient = np.linspace(np.nanmin(data),
                           np.nanmax(data), size)
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
    ax1.set_xlim ((-0.5, size - .5))
    ax1.set_ylim ((-0.5, size - .5))
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
    ax4.vlines(range(size), 0, evect[:,-1], color='k')
    ax4.hlines(0, 0, size, color='red')
    ax4.set_ylabel('E1')
    ax4.set_yticklabels([])
    try:
        ax5.vlines(range(size), 0, evect[:,-2], color='k')
    except IndexError:
        pass
    ax5.hlines(0, 0, size, color='red')
    ax5.set_ylabel('E2')
    ax5.set_yticklabels([])
    try:
        ax6.vlines(range(size), 0, evect[:,-3], color='k')
    except IndexError:
        pass
    ax6.hlines(0, 0, size, color='red')
    ax6.set_ylabel('E3')
    ax6.set_yticklabels([])
    xticklabels = ax4.get_xticklabels() + ax5.get_xticklabels() + ax6.get_xticklabels()
    plt.setp(xticklabels, visible=False)
    if savefig:
        tadbit_savefig(savefig)
    elif show:
        plt.show()
    plt.close('all')


def plot_distance_vs_interactions(data, min_diff=10, max_diff=1000, show=False,
                                  genome_seq=None, resolution=None, axe=None,
                                  savefig=None):
    """
    :param fnam: input file name
    :param 10 min_diff: lower limit (in number of bins)
    :param 1000 max_diff: upper limit (in number of bins) to look for
    :param 100 resolution: group reads that are closer than this resolution
       parameter
    :param None axe: a matplotlib.axes.Axes object to define the plot
       appearance
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    
    """
    resolution = resolution or 1
    dist_intr = dict([(i, 0) for i in xrange(min_diff, max_diff)])
    if isinstance(data, str):
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
                    dist_intr[diff] += 1
                line = fhandler.next()
        except StopIteration:
            pass
        fhandler.close()
    elif isinstance(data, HiC_data):
        max_diff = min(len(data), max_diff)
        if data.section_pos:
            for crm in data.section_pos:
                for diff in xrange(min_diff, min(
                    (max_diff,
                     1 + data.section_pos[crm][1] - data.section_pos[crm][0]))):
                    for i in xrange(data.section_pos[crm][0],
                                    data.section_pos[crm][1] - diff):
                        dist_intr[diff] += data[i, i + diff]
        else:
            for diff in xrange(min_diff, max_diff):
                for i in xrange(len(data) - diff):
                    if not np.isnan(data[i, i + diff]):
                        dist_intr[diff] += data[i,i + diff]
    else:
        if genome_seq:
            max_diff = min(max(genome_seq.values()), max_diff)
            cnt = 0
            for crm in genome_seq:
                for diff in xrange(min_diff, min(
                    (max_diff, genome_seq[crm]))):
                    for i in xrange(cnt, cnt + genome_seq[crm] - diff):
                        if not np.isnan(data[i][i + diff]):
                            dist_intr[diff] += data[i][i + diff]
                cnt += genome_seq[crm]
        else:
            max_diff = min(len(data), max_diff)
            for diff in xrange(min_diff, max_diff):
                for i in xrange(len(data) - diff):
                    if not np.isnan(data[i][i + diff]):
                        dist_intr[diff] += data[i][i + diff]
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
    xp, yp = zip(*sorted(dist_intr.items(), key=lambda x:x[0]))
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
    elif show==True:
        plt.show()


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
    return count_by_len


def insert_sizes(fnam, savefig=None, nreads=None, max_size=99.9, axe=None):
    """
    Plots the distribution of dangling-ends/self-circles lengths
    :param fnam: input file name
    :param None savefig: path where to store the output images.
    :param 99.9 max_size: top percentage of distances to consider, within the
       top 0.01% are usually found very long outliers.
    """
    distr = {}
    genome_seq = OrderedDict()
    fhandler = open(fnam)
    line = fhandler.next()
    while line.startswith('#'):
        if line.startswith('# CRM '):
            crm, clen = line[6:].split()
            genome_seq[crm] = int(clen)
        line = fhandler.next()
    dists = []
    try:
        while True:
            (crm1, pos1, dir1, _, re1, _,
             crm2, pos2, dir2, _, re2) = line.strip().split('\t')[1:12]
            if re1==re2 and crm1 == crm2 and dir1 != dir2:
                dist = abs(int(pos1) - int(pos2))
                dists.append(dist)
                if len(dists) == nreads:
                    break
            line = fhandler.next()
    except StopIteration:
        pass
    fhandler.close()
    ax = setup_plot(axe, figsize=(8, 5))
    max_perc = np.percentile(dists, max_size)
    perc95 = np.percentile(dists, 95)
    perc05 = np.percentile(dists, 5)
    perc99 = np.percentile(dists, 99)
    perc01 = np.percentile(dists, 1)
    ax.hist(dists, bins=100, range=(0, max_perc),
            alpha=0.8, color='darkred')
    ylims = ax.get_ylim()
    plots = []
    plots.append(ax.vlines([perc01], ylims[0], ylims[1],
                      linestyle='--', color='darkgreen',
                      alpha=0.7, linewidth=2))
    plots.append(ax.vlines([perc05], ylims[0], ylims[1],
                      linestyle='--', color='palegreen',
                      alpha=0.7, linewidth=2))
    plots.append(ax.vlines([perc95], ylims[0], ylims[1],
                      linestyle='--', color='orange',
                      alpha=0.7, linewidth=2))
    plots.append(ax.vlines([perc99], ylims[0], ylims[1],
                      linestyle='--', color='orangered',
                      alpha=0.7, linewidth=2))
    ax.set_xlabel('Genomic distance between reads')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of self-circles/dangling-ends ' +
                 'lenghts\n(top %.1f%%, up to %0.f nts)' % (max_size, max_perc))
    plt.subplots_adjust(left=0.1, right=0.75)
    ax.legend(plots, ['1%% (%.0f nts)' % perc01, '5%% (%.0f nts)' % perc05,
                        '95%% (%.0f nts)' % perc95, '99%% (%.0f nts)' % perc99],
              bbox_to_anchor=(1.4, 1), frameon=False)
    if savefig:
        tadbit_savefig(savefig)
    elif not axe:
        plt.show()
    plt.close('all')
    

def plot_genomic_distribution(fnam, first_read=True, resolution=10000,
                              axe=None, ylim=None, savefig=None):
    """
    :param fnam: input file name
    :param True first_read: uses first read.
    :param 100 resolution: group reads that are closer than this resolution
       parameter
    :param None axe: a matplotlib.axes.Axes object to define the plot
       appearance
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    
    """

    distr = {}
    idx1, idx2 = (1, 3) if first_read else (7, 9)
    genome_seq = OrderedDict()
    fhandler = open(fnam)
    line = fhandler.next()
    while line.startswith('#'):
        if line.startswith('# CRM '):
            crm, clen = line[6:].split()
            genome_seq[crm] = int(clen)
        line = fhandler.next()
    try:
        while True:
            crm, pos = line.strip().split('\t')[idx1:idx2]
            pos = int(pos) / resolution
            try:
                distr[crm][pos] += 1
            except KeyError:
                try:
                    distr[crm][pos] = 1
                except KeyError:
                    distr[crm] = {pos: 1}
            line = fhandler.next()
    except StopIteration:
        pass
    fhandler.close()
    if not axe:
        _ = plt.figure(figsize=(15, 3 * len(distr.keys())))

    max_y = max([max(distr[c].values()) for c in distr])
    max_x = max([len(distr[c].values()) for c in distr])
    for i, crm in enumerate(genome_seq if genome_seq else distr):
        plt.subplot(len(distr.keys()), 1, i + 1)
        plt.plot(range(max(distr[crm])),
                 [distr[crm].get(j, 0) for j in xrange(max(distr[crm]))],
                 color='red', lw=1.5, alpha=0.7)
        if ylim:
            plt.vlines(genome_seq[crm] / resolution, ylim[0], ylim[1])
        else:
            plt.vlines(genome_seq[crm] / resolution, 0, max_y)
        plt.xlim((0, max_x))
        plt.ylim(ylim or (0, max_y))
        plt.title(crm)

    if savefig:
        tadbit_savefig(savefig)
    elif not axe:
        plt.show()
    plt.close()


def correlate_matrices(hic_data1, hic_data2, max_dist=10, intra=False,
                       savefig=None, show=False, savedata=None):
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

    :returns: list of correlations and list of genomic distances
    """
    corr = []
    dist = []
    if (intra and hic_data1.sections and hic_data2.sections and 
        hic_data1.sections == hic_data2.sections):
        for i in xrange(1, max_dist + 1):
            diag1 = []
            diag2 = []
            for crm in hic_data1.section_pos:
                for j in xrange(hic_data1.section_pos[crm][0],
                                hic_data1.section_pos[crm][1] - i):
                    diag1.append(hic_data1[j, i + j])
                    diag2.append(hic_data2[j, i + j])
            corr.append(spearmanr(diag1, diag2)[0])
            dist.append(i)
    else:
        if intra:
            warn('WARNING: hic_dta does not contain chromosome coordinates, ' +
                 'intra set to False')
        for i in xrange(1, max_dist + 1):
            diag1 = []
            diag2 = []
            for j in xrange(len(hic_data1) - i):
                diag1.append(hic_data1[j, i + j])
                diag2.append(hic_data2[j, i + j])
            corr.append(spearmanr(diag1, diag2)[0])
            dist.append(i)
    if show or savefig:
        plt.plot(dist, corr, color='orange', linewidth=3, alpha=.8)
        plt.xlabel('Genomic distance in bins')
        plt.ylabel('Spearman rank correlation')
        plt.xlim((0, dist[-1]))
        if savefig:
            tadbit_savefig(savefig)
        if show:
            plt.show()
        plt.close('all')
    if savedata:
        out = open(savedata, 'w')
        out.write('# genomic distance\tSpearman rank correlation\n')
        for i in xrange(len(corr)):
            out.write('%s\t%s\n' % (dist[i], corr[i]))
        out.close()

    return corr, dist


def eig_correlate_matrices(hic_data1, hic_data2, nvect=6,
                           savefig=None, show=False, savedata=None):
    """
    Compare the iteractions of two Hi-C matrices using their 6 first
    eigenvectors, with spearman rank correlation

    :param hic_data1: Hi-C-data object
    :param hic_data2: Hi-C-data object
    :param 6 nvect: number of eigenvectors to compare
    :param None savefig: path to save the plot
    :param False show: displays the plot

    :returns: matrix of correlations
    """
    data1 = hic_data1.get_matrix()
    data2 = hic_data2.get_matrix()
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
        im = axe.imshow(corr, interpolation="nearest",origin='lower')
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

    return corr
