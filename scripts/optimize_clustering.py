"""
17 Jun 2013


"""

from pytadbit import Chromosome
from random import random
from scipy.interpolate import interp1d
import numpy as np
from sys import version_info
if (version_info > (3, 0)):
    from string import ascii_uppercase as uppercase
else:
    from string import uppercase
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
import multiprocessing as mu
from pytadbit.tad_clustering.tad_cmo import optimal_cmo


PATH =  'sample_data/'


def get_distances(tad_matrices, max_num_v=8, n_cpus=8):
    """
    Calculates distances between all pair of tads in the chromosome.
    several CPUs can be used.
    :param tad_matrices: all tads in form {1:[hi-c_list], 2:[hi-c_list]...}
    :param 6 num_v: maximum number of eigenvector to use in the alignment (the
       more the slower... but the better the approximation). Number higher than
       15 should not be considered.
    :param 4 n_cpus: number of CPUs to use
    
    :returns: a dict of distances
    """
    num = len(tad_matrices)
    distances = {}
    cci = {}
    jobs = {}
    pool = mu.Pool(n_cpus)
    for i in xrange(num):
        for j in xrange(i+1, num):
            jobs[(i, j)] = pool.apply_async(
                optimal_cmo, args=(tad_matrices[i], tad_matrices[j]),
                kwds={'max_num_v': max_num_v, 'method': 'frobenius'})
    pool.close()
    pool.join()
    for i in xrange(num):
        for j in xrange(i + 1, num):
            _, _, sc = jobs[(i, j)].get()
            # 0.0001 has shown to be a fair cutoff for p-values for square matrix
            # comparison
            if sc['pval'] < 0.0001:
                cci.setdefault(i, []).append(j)
                distances[(i, j)] = sc['dist']
    return distances, cci


def main():
    """
    main function
    """
    n_pick = 4
    n_tot  = 10
    test_chr = Chromosome(name='Test Chromosome')
    test_chr.add_experiment('exp1', 100000, xp_handler=PATH +
                            'HIC_gm06690_chr19_chr19_100000_obs.txt')
    test_chr.find_tad(['exp1'])
    real_tads = {}
    for i, t in enumerate(test_chr.iter_tads('exp1', normed=False)):
        real_tads[i] = test_chr.experiments['exp1'].tads[i]
        real_tads[i]['hic'] = t[1]
    global DISTRA
    global DISTRD
    DISTRA, DISTRD = get_hic_distr(real_tads)
    # pick some tads
    picked_tads = []
    picked_keys = []
    for i in xrange(n_pick):
        key, new_tad = get_random_tad(real_tads)
        while key in picked_keys or (new_tad['end'] - new_tad['start'] < 15):
            key, new_tad = get_random_tad(real_tads)
        picked_tads.append(new_tad)
        picked_keys.append(key)
    # mutate this tads
    tads = {}
    tad_matrices = []
    tad_names = []
    for i in xrange(n_pick):
        print i
        tads[uppercase[i] + '_' + str(0)] = picked_tads[i]
        tad_names.append(uppercase[i] + '_' + str(0))
        for j in xrange(1, n_tot):
            hic, indels = generate_random_contacts(
                tad1=picked_tads[i]['hic'], prob=0.05, ext=int(random()*4) + 1,
                indel=int(random() * 4) + 1)[1:]
            # indels = '|'.join([str(n-1) if n>0 else '-' + str((abs(n)-1)) for n in indels])
            tads[uppercase[i] + '_' + str(j)] = {
                'hic'  : hic,
                'start': picked_tads[i]['start'],
                'end'  : picked_tads[i]['end']}
            tad_matrices.append(hic)
            tad_names.append(uppercase[i] + '_' + str(j))
    distances, cci = get_distances(tad_matrices, max_num_v=4,
                                   n_cpus=mu.cpu_count())
    results, clusters = pre_cluster(distances, cci, len(tad_matrices))
    paint_clustering(results, clusters, len(tad_matrices), test_chr,
                     tad_names, tad_matrices)


def get_hic_distr(tads):
    # pick all hi-c counts that are not in the diagonal
    all_hic = []
    for t in tads:
        for i, r in enumerate(tads[t]['hic']):
            all_hic += [r[j] for j in range(len(r)) if j != i]
    # all_hic  = [np.log(i) for i in all_hic if i]
    # mean = np.mean(all_hic)
    # std = np.std(all_hic)
    # all_hic  = [(i-mean)/std for i in all_hic]
    all_hic  = [i for i in all_hic if i]
    wina = [0.0]
    cnta = [0.0]
    max_da = np.log(max(all_hic))
    bin_sa = max_da / 100
    binsa = [-bin_sa ] + \
            [float(r)/100 for r in range(0, int(max_da * 100),
                                         int(bin_sa * 100))] + \
            [float('+inf')]
    for b in range(len(binsa)-1):
        wina.append(len([i for i in all_hic if \
                        binsa[b] <= np.log(i) < binsa[b+1] and i]))
        cnta.append(np.exp(binsa[b]))
    try:
        wina = [float(v) / sum(wina) for v in wina]
    except ZeroDivisionError:
        wina = 0.0
    wina = np.cumsum(wina)

    diag_hic = []
    for t in tads:
        for i, r in enumerate(tads[t]['hic']):
            diag_hic += [r[j] for j in range(len(r)) if j == i]
    # diag_hic  = [np.log(i) for i in diag_hic if i]
    # mean = np.mean(diag_hic)
    # std = np.std(diag_hic)
    # diag_hic  = [(i-mean)/std for i in diag_hic]
    diag_hic  = [i for i in diag_hic if i]
    wind = [0.0]
    cntd = [0.0]
    max_dd = np.log(max(diag_hic))
    bin_sd = max_dd / 100
    binsd = [-bin_sd ] + \
            [float(r)/100 for r in range(0, int(max_dd*100),
                                         int(bin_sd*100))] + \
            [float('+inf')]
    for b in range(len(binsd)-1):
        wind.append(len([i for i in diag_hic if \
                        binsd[b] <= np.log(i) < binsd[b+1] and i]))
        cntd.append(np.exp(binsd[b]))
    wind = [float(v) / sum(wind) for v in wind]
    wind = np.cumsum(wind)
    distra = interp1d(wina, cnta, kind='linear')
    distrd = interp1d(wind, cntd, kind='linear')
    return distra, distrd


def generate_1_interaction(diag=False):
    if diag:
        return float(DISTRD(random()))
    if random() > .8: # prob of having a 0.0
        return 0.0
    return float(DISTRA(random()))


def generate_interaction(size):
    rnd_inter = {}
    for i in xrange(size):
        for j in xrange(i, size):
            if i == j:
                rnd_inter[(i, j)] = float(DISTRD(random()))
                continue
            if random() > .8:
                rnd_inter[(i, j)] = 0.1
                rnd_inter[(j, i)] = 0.1
                continue
            rnd_inter[(i, j)] = float(DISTRA(random()))
            rnd_inter[(j, i)] = rnd_inter[(i, j)]
    return rnd_inter


def generate_random_contacts(size=None, tad1=None, prob=0.5,
                             indel=2, ext=5, p_insert=0.5, verbose=False):
    """
    returns 2 random matrices correposnding to TADs
    and 2 contacts lists.
    """
    if not size:
        size1 = size2  = 10# + int(random()*5)
    else:
        size1 = size2  = size# + int(random()*5)
    if tad1:
        size1 = size2 = size = len(tad1)
        contacts1 = {}
        for i in xrange(len(tad1)):
            for j in xrange(len(tad1)):
                contacts1[(i, j)] = tad1[i][j]
    else:
        contacts1 = generate_interaction(size1)
    contacts2 = deepcopy(contacts1)
    # randomly change some contacts
    for i in xrange(size1):
        for j in xrange(i, size1):
            if random() < prob:
                contacts2[(i, j)] = generate_1_interaction(i==j)
                contacts2[(j, i)] = contacts2[(i, j)]
    # a bit removing and inserting columns
    changes = 0
    p_delete = 1 - p_insert
    indels = []
    while changes < indel:
        rnd = random()
        my_j = int(random()*size2)
        if rnd < p_insert:
            for _ in xrange(ext):
                indels.append(my_j + 1)
                size2 += 1
                for i, d in enumerate(indels):
                    if abs(d) > (my_j + 1):
                        indels[i] = (d + 1) if d > 0 else -(abs(d) + 1)
                for i in xrange(size2 - 2, -1, -1):
                    for j in xrange(size2 - 2, i - 1, -1):
                        if i >= my_j and j >= my_j:
                            contacts2[(i + 1, j + 1)] = contacts2[(j, i)]
                        if j >= my_j:
                            contacts2[(i, j + 1)]     = contacts2[(j, i)]
                    if i == my_j:
                        contacts2[(i, j + 1)]     = contacts2[(j, i)]
                for i in xrange(size2):
                    for j in xrange(i, size2):
                        contacts2[(j, i)] = contacts2[(i, j)]
                for i in xrange(size2):
                    contacts2[(i, my_j)] = generate_1_interaction(i==my_j)
                    contacts2[(my_j, i)] = contacts2[(i, my_j)]
        elif rnd > 1 - p_delete and size2 > ext:
            my_j = min(my_j, size2 - ext)
            for _ in xrange(ext):
                size2 -= 1
                indels.append(-(my_j + 1))
                for i, d in enumerate(indels):
                    if abs(d) > (my_j + 1):
                        indels[i] = d - 1  if d > 0 else -(abs(d) - 1)
                for i in xrange(size2):
                    for j in xrange(i, size2):
                        if i >= my_j and i >= my_j:
                            contacts2[(i, j)] = contacts2[(i + 1, j + 1)]
                        elif j >= my_j:
                            contacts2[(i, j)] = contacts2[(i , j + 1)]
                for i in xrange(size2):
                    for j in xrange(i + 1, size2):
                        contacts2[(j, i)] = contacts2[(i, j)]
                for i in xrange(size2):
                    del(contacts2[(i, size2)])
                    del(contacts2[(size2, i)])
                del(contacts2[(size2, size2)])
        else:
            continue
        changes += 1
    print 'inserted/deleted', [str(i-1) if i>0 else '-' + str((abs(i)-1)) \
                               for i in indels]
    for i in range(size1):
        if not (i, i) in contacts1:
            contacts1[(i, i)] = generate_1_interaction(True)
    for i in range(size2):
        if not (i, i) in contacts2:
            contacts2[(i, i)] = generate_1_interaction(True)
    tad1 = contact2matrix(contacts1, size1)
    tad2 = contact2matrix(contacts2, size2)
    if verbose:
        print '\t'.join(['\033[0;31m{}\033[0m'.format(i) for i in range(len(tad1))])# + '\n'
        for t in tad1: print '\t'.join(['{:.1f}'.format(i) for i in t])# + '\n'
        print '-'*20
        print '\t'.join(['\033[0;31m{}\033[0m'.format(i) for i in range(len(tad2))])# + '\n'
        for t in tad2: print '\t'.join(['{:.1f}'.format(i) for i in t])# + '\n'
    return tad1, tad2, indels#, contacts1, contacts2


def paint_clustering(results, clusters, num, chrom, tad_names, matrix):
    dendros = []
    axes = []
    prev = 0
    xlim = [-100, 100]
    tmp = []
    for i, result in enumerate(results):
        if axes:
            axes[-1].set_xticklabels([], visible=False)
        clust = linkage(result, method='ward')
        tmp = dendrogram(clust, orientation='right', no_plot=True)['leaves']
        dendros += reversed(list([clusters[i][n] for n in tmp]))
        axes.append(plt.subplot2grid((num, 9),(prev, 0), rowspan=len(result),
                                     colspan=4))
        dendrogram(clust, orientation='right',
                   labels=[tad_names[c] for c in clusters[i]])
        if xlim[0] < axes[-1].get_xlim()[0]:
            xlim[0] = axes[-1].get_xlim()[0]
        if xlim[1] > axes[-1].get_xlim()[1]:
            xlim[1] = axes[-1].get_xlim()[1]
        prev += len(result)
    for ax in axes:
        ax.set_xlim(left=xlim[0], right=xlim[1])
    axes = []
    for i, j in enumerate(dendros):
        axes.append(plt.subplot2grid((num, 9),(i, 4)))#gs1[i]))
        chrom.visualize('exp1',
                        tad=matrix[j],
                        axe=axes[-1], show=False)
        axes[-1].set_axis_off()
    ax4 = plt.subplot2grid((num, 9),(0, 5), rowspan=num, colspan=4)
    chrom.visualize('exp1', paint_tads=True, axe=ax4)
    plt.draw()


def contact2matrix(contacts, size):
    matrix = [[0  for _ in xrange(size)] for _ in xrange(size)]
    for i, j in contacts:
        try:
            matrix[i][j] = contacts[(i, j)]
            matrix[j][i] = contacts[(j, i)]
        except IndexError:
            pass
    return matrix


def get_random_tad(tads):
    tad = tads.keys()[int(len(tads) * random())]
    return tad, tads[tad]


def pre_cluster(distances, cci, num):
    """
    """
    clusters = []
    list_nodes = cci.keys()
    while list_nodes:
        node = list_nodes.pop()
        found = False
        for cluster in clusters:
            if node in cluster:
                found = True
                cluster = list(set(cluster+cci.get(node, [])))
        if not found:
            clusters.append(list(set([node] + cci.get(node, []))))
        while True:
            for i, cluster1 in enumerate(clusters):
                for j, cluster2 in enumerate(clusters[i+1:]):
                    if set(cluster1) & set(cluster2):
                        clusters[i] = list(set(cluster1+cluster2))
                        del(clusters[i+j+1])
                        break
                else:
                    continue
                break
            else:
                break
    results = [[[0 for _ in xrange(len(i))] for _ in xrange(len(i))] \
               for i in clusters]
    for k, cluster in enumerate(clusters):
        trans = [i for i in xrange(num) if i in cluster]
        for i in xrange(num):
            for j in xrange(i + 1, num):
                if not (i in cluster and j in cluster):
                    continue
                if not (i, j) in distances:
                    continue
                results[k][trans.index(i)][trans.index(j)] = distances[(i, j)]
                results[k][trans.index(j)][trans.index(i)] = distances[(i, j)]
    return results, clusters


if __name__ == "__main__":
    exit(main())
