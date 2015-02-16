"""
18 Nov 2014


"""
from pytadbit.utils.extraviews import tadbit_savefig
from pytadbit.utils.tadmaths import nozero_log_matrix as nozero_log
from warnings import warn
from collections import OrderedDict
import numpy as np
from pytadbit.parsers.hic_parser import load_hic_data_from_reads
from scipy.stats import norm as sc_norm, skew, kurtosis
from scipy.stats import pearsonr, spearmanr
from numpy.linalg import eigh
import os

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')

def hic_map(data, resolution=None, normalized=False, masked=None,
            by_chrom=False, savefig=None, show=False, savedata=None,
            focus=None, clim=None, cmap='Reds', pdf=False, decay=False):
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
    :param False decay: plot the correlation between genomic distance and
       interactions (usually a decay).
    :param None clim: cutoff for the upper and lower bound in the coloring scale
       of the heatmap
    :param False pdf: when using the bny_chrom option, to specify the format of
       the stored images
    :param Reds cmap: color map to be used for the heatmap
    """
    if isinstance(data, str):
        data = load_hic_data_from_reads(data, resolution=resolution)
    hic_data = data
    if hic_data.bads and not masked:
        masked = hic_data.bads
    # save and draw the data
    if by_chrom:
        if focus:
            raise Exception('Incompatible options focus and by_chrom\n')
        os.system('mkdir -p ' + savedata)
        if not clim:
            clim = (np.log2(min(hic_data.values())),
                    np.log2(max(hic_data.values())))
        for i, crm1 in enumerate(hic_data.chromosomes):
            for crm2 in hic_data.chromosomes.keys()[i:]:
                if by_chrom == 'intra' and crm1 != crm2:
                    continue
                if by_chrom == 'inter' and crm1 == crm2:
                    continue
                subdata = hic_data.get_matrix(focus=(crm1, crm2), normalized=normalized)
                if savedata:
                    out = open('%s/%s.mat' % (
                        savedata, '_'.join(set((crm1, crm2)))), 'w')
                    out.write('\n'.join(['\t'.join([str(i) for i in d])
                                         for d in subdata]) + '\n')
                    out.close()
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
                             resolution=resolution)
    else:
        if savedata:
            out = open(savedata, 'w')
            out.write('\n'.join(
                ['\t'.join([str(i) for i in line])
                 for line in hic_data.get_matrix(
                     focus=focus, normalized=normalized)]) + '\n')
            out.close()
        if show or savefig:
            subdata = hic_data.get_matrix(focus=focus, normalized=normalized)
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
                     clim=clim, cmap=cmap, resolution=resolution)


def draw_map(data, genome_seq, cumcs, savefig, show, resolution=None, one=False,
             clim=None, cmap='Reds', decay=False):
    _ = plt.figure(figsize=(15.,12.5))
    ax1 = plt.axes([0.34, 0.08, 0.6, 0.7205])
    ax2 = plt.axes([0.07, 0.65, 0.21, 0.15])
    if decay:
        ax3 = plt.axes([0.07, 0.42, 0.21, 0.15])
        plot_distance_vs_interactions(data, axe=ax3, resolution=resolution)
    ax4 = plt.axes([0.34, 0.805, 0.6, 0.04], sharex=ax1)
    ax5 = plt.axes([0.34, 0.845, 0.6, 0.04], sharex=ax1)
    ax6 = plt.axes([0.34, 0.885, 0.6, 0.04], sharex=ax1)
    cmap = plt.get_cmap(cmap)
    cmap.set_bad('darkgrey', 1)
    data = nozero_log(data, np.log2)
    ax1.imshow(data, interpolation='none',
               cmap=cmap, vmin=clim[0] if clim else None,
               vmax=clim[1] if clim else None)
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
    h  = ax2.hist(data, color='black', linewidth=2,
                   bins=20, histtype='step', normed=True)
    _  = ax2.imshow(gradient, aspect='auto', cmap=cmap,
                     extent=(np.nanmin(data), np.nanmax(data) , 0, max(h[0])))
    if genome_seq:
        for crm in genome_seq:
            ax1.vlines(cumcs[crm][0]-.5, cumcs[crm][0]-.5, cumcs[crm][1]-.5, color='k',
                       linestyle=':')
            ax1.vlines(cumcs[crm][1]-.5, cumcs[crm][0]-.5, cumcs[crm][1]-.5, color='k',
                       linestyle=':')
            ax1.hlines(cumcs[crm][1]-.5, cumcs[crm][0]-.5, cumcs[crm][1]-.5, color='k',
                       linestyle=':')
            ax1.hlines(cumcs[crm][0]-.5, cumcs[crm][0]-.5, cumcs[crm][1]-.5, color='k',
                       linestyle=':')

        if not one:
            vals = [0]
            keys = ['']
            for crm in genome_seq:
                vals.append(cumcs[crm][0])
                keys.append(crm)
            vals.append(cumcs[crm][1])
            ax1.set_yticks(vals)
            ax1.set_yticklabels('')
            ax1.set_yticks([float(vals[i]+vals[i+1])/2 for i in xrange(len(vals) - 1)], minor=True)
            ax1.set_yticklabels(keys, minor=True)
            for t in ax1.yaxis.get_minor_ticks():
                t.tick1On = False
                t.tick2On = False 
    ax2.set_xlim((np.nanmin(data), np.nanmax(data)))
    ax2.set_ylim((0, max(h[0])))
    ax1.set_xlim ((-0.5, size - .5))
    ax1.set_ylim ((-0.5, size - .5))
    ax2.set_xlabel('log interaction count')
    normfit = sc_norm.pdf(data, np.mean(data), np.std(data))
    ax2.plot(data, normfit, 'w.', markersize=2.5, alpha=.4)
    ax2.plot(data, normfit, 'k.', markersize=1.5, alpha=1)
    ax2.set_title('skew: %.3f, kurtosis: %.3f' % (skew(data),
                                                   kurtosis(data)))
    ax4.vlines(range(size), 0, evect[:,-1], color='k')
    ax4.hlines(0, 0, size, color='red')
    ax4.set_ylabel('E1')
    ax4.set_yticklabels([])
    plt.setp(ax4, 'xticklabels', [])
    ax5.vlines(range(size), 0, evect[:,-2], color='k')
    ax5.hlines(0, 0, size, color='red')
    ax5.set_ylabel('E2')
    ax5.set_yticklabels([])
    plt.setp(ax5, 'xticklabels', [])
    ax6.vlines(range(size), 0, evect[:,-3], color='k')
    ax6.hlines(0, 0, size, color='red')
    ax6.set_ylabel('E3')
    ax6.set_yticklabels([])
    plt.setp(ax6, 'xticklabels', [])
    if savefig:
        tadbit_savefig(savefig)
    elif show:
        plt.show()
    plt.close('all')


def plot_distance_vs_interactions(data, min_diff=10, max_diff=1000000,
                                  resolution=None, axe=None, savefig=None):
    """
    :param fnam: input file name
    :param 100 min_diff: lower limit kn genomic distance (usually equal to read
       length)
    :param 1000000 max_diff: upper limit in genomic distance to look for
    :param 100 resolution: group reads that are closer than this resolution
       parameter
    :param None axe: a matplotlib.axes.Axes object to define the plot
       appearance
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    
    """
    resolution = resolution or 1
    dist_intr = dict([(i, 0.00001) for i in xrange(max_diff)])
    if isinstance(data, str):
        fhandler = open(data)
        line = fhandler.next()
        while line.startswith('#'):
            line = fhandler.next()
        try:
            while True:
                _, cr1, ps1, _, _, _, _, cr2, ps2, _ = line.rsplit('\t', 9)
                if cr1 != cr2:
                    line = fhandler.next()
                    continue
                diff = abs(int(ps1) - int(ps2))
                if max_diff > diff > min_diff:
                    dist_intr[diff] += 1
                line = fhandler.next()
        except StopIteration:
            pass
        fhandler.close()
        for k in dist_intr.keys()[:]:
            if dist_intr[k] <= 2:
                del(dist_intr[k])
    else:
        max_diff = min(len(data), max_diff)
        dist_intr = dict([(i, 0.00001) for i in xrange(1, max_diff)])
        for diff in xrange(1, max_diff):
            for i in xrange(len(data) - diff):
                if not np.isnan(data[i][i + diff]):
                    dist_intr[diff] += np.exp2(data[i][i + diff])
    if not axe:
        fig=plt.figure()
        axe = fig.add_subplot(111)
    x, y = zip(*sorted(dist_intr.items(), key=lambda x:x[0]))
    axe.plot(x, y, 'k.')
    coeffs = np.polyfit(np.log(x[:5*len(x)/10]),
                        np.log(y[:5*len(x)/10]), 1)
    poly = np.poly1d(coeffs)
    yfit = lambda x: np.exp(poly(np.log(x)))
    # plot line of best fit
    axe.plot(x,yfit(x), color= 'red', lw=2)
    axe.set_ylabel('Log interaction count')
    axe.set_xlabel('Log genomic distance (binned by %d bp)' % resolution)
    axe.set_xlim((np.log(min_diff), np.log(max_diff)))
    axe.set_title('Slope (first half, up to %d bins)\n' % (max_diff / 2) +
                  r'$\alpha=%.2f$' % (coeffs[0]))
    axe.set_xscale('log')
    axe.set_yscale('log')
    axe.set_xlim((0, max_diff))
    # plt.legend()
    if savefig:
        tadbit_savefig(savefig)
    elif not axe:
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
    for i, fnam in enumerate([fnam1, fnam2]):
        fhandler = open(fnam)
        line = fhandler.next()
        while line.startswith('#'):
            line = fhandler.next()
        try:
            count_by_len[i] = {}
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
    plt.xlabel('read length (bp)')
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


def plot_genomic_distribution(fnam, first_read=True, resolution=10000,
                              axe=None, ylim=None, savefig=None):
    """
    :param fnam: input file name
    :param True first_read: map first read.
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


def correlate_matrices(hic_data1, hic_data2, resolution=1, max_dist=10,
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
    :param False show: displays the plot

    :returns: list of correlations and list of genomic distances
    """
    corr = []
    dist = []
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
    corr = []
    ev1, evect1 = eigh(np.array([[hic_data1[i, j]
                         for j in xrange(len(hic_data1))]
                        for i in xrange(len(hic_data1))]))
    ev2, evect2 = eigh(np.array([[hic_data2[i, j]
                         for j in xrange(len(hic_data2))]
                        for i in xrange(len(hic_data2))]))
    corr = [[0 for _ in xrange(nvect)] for _ in xrange(nvect)]
    sort_perm = ev1.argsort()
    evect1 = evect1[sort_perm][::-1]
    sort_perm = ev2.argsort()
    evect2 = evect2[sort_perm][::-1]
    for i in xrange(nvect):
        for j in xrange(nvect):
            corr[i][j] = abs(pearsonr(evect1[:,i],
                                      evect2[:,j])[0])
    if show or savefig:
        plt.imshow(corr, interpolation="nearest",origin='lower')
        plt.xlabel('Eigen Vectors exp. 1')
        plt.ylabel('Eigen Vectors exp. 2')
        plt.xticks(range(nvect), range(1, nvect + 1))
        plt.yticks(range(nvect), range(1, nvect + 1))
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Pearson correlation')
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
