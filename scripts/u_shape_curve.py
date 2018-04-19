"""
14 Jun 2013

Sample script to infer bias toward TAD boundaries

"""

from pytadbit import Chromosome
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from matplotlib import pyplot as plt

PATH =  'sample_data/'


def tad_breaker(tads, cluster, resolution, bins=20, show_plot=False,
                title=None, max_tad_size=3000000):
    """
    Find out if TAD boundaries are breaking clusters of colocalizing
    (epi)genetic elements.

    :param tads:
    :param cluster: cluster of elements
    :param title: for plotting
    :param resolution: in bases
    :param 20 bins: number of bins to use
    :param False show_plot:
    :para False test: for significant difference between extremes and center

    :returns: list of counts for each bin
    """
    def arrange(beg, end, num):
        num -= 1
        step  = int((end - beg) / num)
        for k in xrange(int(beg), int(end), step):
            yield k
        yield int(end)
    brk_set = []
    for t in xrange(len(tads)):
        if not tads[t]['score'] >= 5:
           continue
        if (tads[t]['end'] - tads[t]['start']) * resolution > max_tad_size:
            continue
        brk_set.append([b for b in arrange(tads[t]['start'] * resolution,
                                           tads[t]['end'] * resolution, bins)])
    counts = []
    for b in range(bins):
        count = 0
        for start, end in cluster:
            for brk in brk_set:
                if start < brk[b] < end:
                    count += 1
                    break
        counts.append(1-float(count)/len(cluster))
    print counts
    if show_plot:
        plt.bar(range(bins), counts)
        plt.title(title)
        plt.ylim((0.8, 1))
        plt.show()
    return counts


def get_genes():
    """
    list of genes annotated with at least one of these GO term:
    GO:0043021
    GO:0043228
    GO:0031974
    GO:0003723
    GO:0006396
    GO:0003735
    GO:0005840
    GO:0005739
    GO:0007049
    (dixon paper)
    """
    crms = {}
    geneids = {}
    for line in open(PATH + 'human_cellcomp_ensembl_v54.tsv'):
        if not line.startswith('ENS'):
            continue
        geneid, crm, beg, _ = line.split()
        if len(crm) > 2:
            continue
        geneids.setdefault(crm, []).append([geneid, float(beg)])
        crms.setdefault(crm, []).append(float(beg))
    distmatrix = {}
    for crm in crms:
        distmatrix[crm] = [[0 for _ in xrange(len(crms[crm]))]
                           for _ in xrange(len(crms[crm]))]
        for i in xrange(len(crms[crm])):
            for j in xrange(len(crms[crm])):
                distmatrix[crm][i][j] = abs(crms[crm][i] - crms[crm][j])
    return distmatrix, geneids


def main():
    """
    main function
    """
    # retieve HOX genes

    distmatrix, geneids = get_genes()
    # compute TADs for human chromosome 19
    test_chr = Chromosome(name='Test Chromosome')
    test_chr.add_experiment('exp1', 100000, xp_handler=PATH +
                            'HIC_gm06690_chr19_chr19_100000_obs.txt')
    test_chr.find_tad(['exp1'])
    exp = test_chr.experiments['exp1']
    clust = linkage(distmatrix['19'])
    cl_idx = list(fcluster(clust, t=1, criterion='inconsistent'))
    print max(cl_idx), 'clusters'
    cluster=[[] for _ in xrange(1, max(cl_idx)+1)]
    for i, j in enumerate(cl_idx):
        cluster[j-1].append(geneids['19'][i][1])
    for i, _ in enumerate(cluster):
        cluster[i] = min(cluster[i]), max(cluster[i])
    tad_breaker(exp.tads, cluster, exp.resolution, show_plot=True, bins=5,
                title='Proportion of HOX genes according to position in a TAD')



if __name__ == "__main__":
    exit(main())
