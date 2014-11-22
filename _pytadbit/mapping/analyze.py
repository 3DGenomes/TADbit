"""
18 Nov 2014


"""
from pytadbit.utils.extraviews import tadbit_savefig
from warnings import warn
import numpy as np
# from scipy.ndimage import gaussian_filter1d

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')

def hic_map(fnam, genome_seq, masked=None, resolution=100000,
            savefig=None):
    cumcs = {} 
    total = 0
    for crm in genome_seq:
        cumcs[crm] = total
        total += len(genome_seq[crm]) / resolution
    # bin the data
    data = [[0 for _ in xrange(total + 1)] for _ in xrange(total + 1)]
    masked = masked or set()
    for line in open(fnam):
        read, cr1, ps1, _, _, _, _, cr2, ps2, _, _, _, _ = line.split()
        if read in masked:
            continue
        ps1 = int(ps1) / resolution
        ps2 = int(ps2) / resolution
        try:
            data[cumcs[cr1] + ps1][cumcs[cr2] + ps2] += 1
        except:
            break
    # do the plot
    import numpy as np
    plt.figure(figsize=(16, 12))
    plt.imshow(np.log2(data), origin='lower', cmap='gist_earth')
    plt.colorbar()
    if savefig:
        tadbit_savefig(savefig)
    else:
        plt.show()

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
    for line in open(fnam):
        _, cr1, ps1, _, _, _, _, cr2, ps2, _ = line.rsplit('\t', 9)
        if cr1 != cr2:
            continue
        diff = resolution * (abs(int(ps1) - int(ps2)) / resolution)
        if max_diff > diff > min_diff:
            dist_intr.setdefault(diff, 0)
            dist_intr[diff] += 1
            
    for k in dist_intr.keys()[:]:
        if dist_intr[k] <= 2:
            del(dist_intr[k])
                    
    if not axe:
        fig=plt.figure()
        ax = fig.add_subplot(111)

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
        ax = fig.add_subplot(111)
    colors = ['olive', 'darkcyan']
    for i, fnam in enumerate([fnam1, fnam2]):
        print fnam
        count_by_len[i] = {}
        for line in open(fnam):
            _, length, _, _ = line.rsplit('\t', 3)
            try:
                count_by_len[i][int(length)] += 1
            except KeyError:
                count_by_len[i][int(length)] = 1
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
    for line in open(fnam):
        crm, pos = line.split()[idx1:idx2]
        pos = int(pos) / resolution
        try:
            distr[crm][pos] += 1
        except KeyError:
            try:
                distr[crm][pos] = 1
            except KeyError:
                distr[crm] = {pos: 1}

    if not axe:
        fig=plt.figure(figsize=(15, 3 * len(distr.keys())))

    max_y = max([max(distr[c].values()) for c in distr])
    max_x = max([len(distr[c].values()) for c in distr])
    for i, crm in enumerate(distr):
        plt.subplot(len(distr.keys()), 1, i + 1)
        plt.plot(range(max(distr[crm])),
                 [distr[crm].get(j, 0) for j in xrange(max(distr[crm]))],
                 color='red', lw=1.5, alpha=0.7)
        plt.xlim((0, max_x))
        plt.ylim(ylim or (0, max_y))
        plt.title(crm)

    if savefig:
        tadbit_savefig(savefig)
    elif not axe:
        plt.show()

