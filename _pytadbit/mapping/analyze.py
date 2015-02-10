"""
18 Nov 2014


"""
from pytadbit.utils.extraviews import tadbit_savefig
from warnings import warn
from collections import OrderedDict
import numpy as np
from scipy.stats import norm as sc_norm, skew, kurtosis
import os

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')

def hic_map(data, biases=None, masked=None, resolution=100000, by_chrom=False,
            savefig=None, show=False, savedata=None, focus=None, cmap='Reds'):
    if isinstance(data, str):
        genome_seq = OrderedDict()
        fhandler = open(data)
        line = fhandler.next()
        while line.startswith('#'):
            if line.startswith('# CRM '):
                crm, clen = line[6:].split()
                genome_seq[crm] = int(clen)
            line = fhandler.next()
        cumcs = {} 
        total = 0
        for crm in genome_seq:
            cumcs[crm] = (total, total + genome_seq[crm] / resolution + 1)
            total += genome_seq[crm] / resolution + 1
        # bin the data
        data = [[0 for _ in xrange(total + 1)] for _ in xrange(total + 1)]
        masked = masked or set()
        try:
            while True:
                read, cr1, ps1, _, _, _, _, cr2, ps2, _ = line.split('\t', 9)
                if read in masked:
                    continue
                ps1 = int(ps1) / resolution
                ps2 = int(ps2) / resolution
                try:
                    data[cumcs[cr1][0] + ps1][cumcs[cr2][0] + ps2] += 1
                    data[cumcs[cr2][0] + ps2][cumcs[cr1][0] + ps1] += 1
                except:
                    break
                line = fhandler.next()
        except StopIteration:
            pass
        fhandler.close()
    else:
        hic_data = data
        beg, end = focus if focus else (0, len(hic_data))
        beg -= 1 if focus else 0
        if biases:
            data = [[hic_data[len(hic_data) * i + j] / (biases[i] * biases[j])
                     for j in xrange(beg, end)]
                    for i in xrange(beg, end)]
        else: 
            data = [[hic_data[len(hic_data) * i + j]
                     for j in xrange(beg, end)]
                    for i in xrange(beg, end)]
        if masked:
            data = [[float('nan') for j in xrange(len(hic_data))]
                    if i in masked else data[i] for i in xrange(len(hic_data))]
            data = [[float('nan') if j in masked else data[i][j]
                     for j in xrange(len(hic_data))]
                    for i in xrange(len(hic_data))]
    # save the data
    if savedata:
        if by_chrom:
            if focus:
                raise Exception('Incompatible options focus and by_chrom\n')
            os.system('mkdir -p ' + savedata)
            for i, crm1 in enumerate(genome_seq):
                for crm2 in genome_seq.keys()[i:]:
                    if by_chrom == 'intra' and crm1 != crm2:
                        continue
                    if by_chrom == 'inter' and crm1 == crm2:
                        continue
                    out = open('%s/%s.mat' % (savedata,
                                              '_'.join(set((crm1, crm2)))), 'w')
                    subdata = [[cell for cell in
                                line[cumcs[crm2][0]:cumcs[crm2][1]]]
                               for line in data[cumcs[crm1][0]:cumcs[crm1][1]]]
                    out.write('\n'.join(['\t'.join([str(i) for i in d])
                                         for d in subdata]) + '\n')
                    out.close()
                    if show or savefig:
                        draw_map(subdata, 
                                 OrderedDict([(k, genome_seq[k])
                                              for k in genome_seq.keys()[:]
                                              if k in [crm1, crm2]]),
                                 resolution,
                                 '%s/%s.pdf' % (savefig,
                                                '_'.join(set((crm1, crm2)))),
                                 show, cmap=cmap)
        else:
            out = open(savedata, 'w')
            for line in data:
                out.write('\t'.join([str(cell) for cell in line]) + '\n')
            out.close()
            if show or savefig:
                draw_map(data, genome_seq, resolution, savefig, show, cmap='Reds')

def draw_map(data, genome_seq, resolution, savefig, show, cmap='Reds'):
    cumcs = {} 
    total = 0
    for crm in genome_seq:
        cumcs[crm] = (total, total + genome_seq[crm] / resolution + 1)
        total += genome_seq[crm] / resolution + 1
    fig = plt.figure(figsize=(10.,9.1))
    axe = fig.add_subplot(111)
    axe.imshow(np.log2(data), interpolation='none',
               cmap=cmap)
    fig.subplots_adjust(top=.8, left=0.35)
    size = len(data)
    data = np.log2([(i or 0.1)  for j in data for i in j])
    gradient = np.linspace(min(data), max(data), size)
    gradient = np.vstack((gradient, gradient))
    axe2 = fig.add_axes([0.1, 0.78, 0.2, 0.1])
    h  = axe2.hist(data, color='grey', bins=20, histtype='step', normed=True)
    _  = axe2.imshow(gradient, aspect='auto', cmap=cmap,
                     extent=(min(data), max(data) , 0, max(h[0])))
    for crm in genome_seq:
        axe.vlines(cumcs[crm][0], cumcs[crm][0], cumcs[crm][1], color='k',
                   linestyle=':')
        axe.vlines(cumcs[crm][1], cumcs[crm][0], cumcs[crm][1], color='k',
                   linestyle=':')
        axe.hlines(cumcs[crm][1], cumcs[crm][0], cumcs[crm][1], color='k',
                   linestyle=':')
        axe.hlines(cumcs[crm][0], cumcs[crm][0], cumcs[crm][1], color='k',
                   linestyle=':')
    axe2.set_xlim((min(data), max(data)))
    axe2.set_ylim((0, max(h[0])))
    axe.set_xlim ((-0.5, size - .5))
    axe.set_ylim ((-0.5, size - .5))
    axe2.set_xlabel('log interaction count', size='small')
    normfit = sc_norm.pdf(data, np.mean(data), np.std(data))
    axe2.plot(data, normfit, 'g.', markersize=1, alpha=.1)
    axe2.set_title('skew: %.3f, kurtosis: %.3f' % (skew(data),
                                                   kurtosis(data)),
                   size='small')
    vals = [0]
    keys = ['']
    for crm in genome_seq:
        vals.append(cumcs[crm][0])
        keys.append(crm)
    vals.append(cumcs[crm][1])
    axe.set_yticks(vals)
    axe.set_yticklabels('')
    axe.set_yticks([float(vals[i]+vals[i+1])/2 for i in xrange(len(vals) - 1)], minor=True)
    axe.set_yticklabels(keys, minor=True)
    for t in axe.yaxis.get_minor_ticks():
        t.tick1On = False
        t.tick2On = False 
    if savefig:
        tadbit_savefig(savefig)
    elif show:
        plt.show()
    plt.close('all')

def plot_distance_vs_interactions(fnam, min_diff=100, max_diff=1000000,
                                  resolution=100, axe=None, savefig=None):
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
    dist_intr = {}
    fhandler = open(fnam)
    line = fhandler.next()
    while line.startswith('#'):
        line = fhandler.next()
    try:
        while True:
            _, cr1, ps1, _, _, _, _, cr2, ps2, _ = line.rsplit('\t', 9)
            if cr1 != cr2:
                line = fhandler.next()
                continue
            diff = resolution * (abs(int(ps1) - int(ps2)) / resolution)
            if max_diff > diff > min_diff:
                dist_intr.setdefault(diff, 0)
                dist_intr[diff] += 1
            line = fhandler.next()
    except StopIteration:
        pass
    fhandler.close()
            
    for k in dist_intr.keys()[:]:
        if dist_intr[k] <= 2:
            del(dist_intr[k])
                    
    if not axe:
        fig=plt.figure()
        _ = fig.add_subplot(111)
        
    x, y = zip(*sorted(dist_intr.items(), key=lambda x:x[0]))
    plt.plot(x, y, 'k.')
    # sigma = 10
    # p_x = gaussian_filter1d(x, sigma)
    # p_y = gaussian_filter1d(y, sigma)
    # plot line of best fit
    # plt.plot(p_x, p_y,color= 'darkgreen', lw=2, label='Gaussian fit')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Log genomic distance (binned by %d bp)' % resolution)
    plt.ylabel('Log interaction count')
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

